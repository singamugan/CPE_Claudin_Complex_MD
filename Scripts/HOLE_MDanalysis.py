#!/usr/bin/env python
# coding: utf-8

import sys
import MDAnalysis as mda


from MDAnalysis.analysis import hole2
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os
from ast import literal_eval
import warnings

###############################################################################################################################################



def hole_analysis(u, u_last_frame, frames, channel_ax, new_path):


    ## The MDAnalysis.analysis.hole2 module provides wrapper classes and functions that call the HOLE program.
    ## This means you must have installed the program yourself before you can use the class.
    ## You then need to pass the path to your hole executable to the class; in this script, hole is installed at ~/hole2/exe/hole.


    print("Beginning HOLE analysis of the full trajectory...\n")
    
    
    if channel_ax == "x":
    	channel_vec = [1,0,0]
    elif channel_ax == "y":
    	channel_vec = [0,1,0]
    elif channel_ax == "z":
    	channel_vec = [0,0,1]
    
    
    ha = hole2.HoleAnalysis(u, select='protein', cvect=channel_vec, executable='~/hole2/exe/hole')
    ha.run()
   
    ## HOLE writes out files for each frame.   
    ## The data is stored in ha.profiles as a dictionary of numpy.recarrays .
    ## The dictionary is indexed by frame.
    
    
    
    ## plotting the minimum pore radius over all profiles, i.e for each frame
    
    min_radii = ha.min_radius()
    
    print("Saving plot of the minimum pore radius for each frame to PNG-file and to CSV-file.\n")
    
    plt.plot(min_radii[:, 0], min_radii[:, 1])
    plt.ylabel('Minimum HOLE radius $R$ ($\AA$)')
    plt.xlabel('Frame')
    plt.title("Minimum pore radius of all frames")
    plt.savefig(new_path+'min_pore_radius_over_all_frames.png')   ## save the figure to file
    plt.close()    ## close the figure window
    
    
    data_min = {'Frame': min_radii[:, 0], 'Minimum HOLE radius': min_radii[:, 1]}
    df_min = pd.DataFrame(data=data_min)
    df_min.to_csv(new_path+'min_pore_radius_over_all_frames.csv', index=False)
    
    
    ###################################################################################################
    
    ## plotting the mean of the pore radius over the slected frame range by
    ## collecting the pore radii into bins by reaction
    
    frame_range = range((frames[0]-1),frames[1],1)
    n_std=1	# how many standard deviations from mean
    
    binned, bins = ha.bin_radii(frames=frame_range, bins=100, range=None)
    mean = np.array(list(map(np.mean, binned)))
    midpoints = 0.5 * (bins[1:] + bins[:-1])
    
    print("Saving plot of the mean of the pore radius over the pore coordinate to PNG-file and to CSV-file for the selected frame range.\n")
    
    fig, ax = plt.subplots()
    std = np.array(list(map(np.std, binned)))
    ax.fill_between(midpoints, mean-(n_std*std), mean+(n_std*std), color='blue', alpha=0.2, label='{} std'.format(n_std))
    ax.plot(midpoints, mean, color='blue', linestyle='-', label='mean')
    ax.set_xlabel(r"Pore coordinate $\zeta$ ($\AA$)")
    ax.set_ylabel(r"HOLE radius $R$ ($\AA$)")
    ax.legend(loc='best')
    plt.title("Mean and standard deviation of the pore radius \nover the pore coordinate for frame "+ str(frames[0]) + '-' + str(frames[1]))
    
    fig.savefig(new_path+"mean_pore_radius_over_pore_coord.png")   ## save the figure to file
    plt.close(fig)    ## close the figure window
    
    data_mean = {'Pore coordinate': midpoints, 'Mean HOLE radius': mean}
    df_mean = pd.DataFrame(data=data_mean)
    df_mean['S.D.'] = std
    df_mean.to_csv('mean_pore_radius_over_pore_coord.csv', index=False)
    
    ###################################################################################################

 
    ## To write pore surfaces to a VMD-file, we can use the 'create_vmd_surface()' which is built into the HoleAnalysis class.
    ## However, for universes with a huge number of frames this might not work properly,
    ## which is why we simply write a VMD file only for the last frame.
    
    print("Writing VMD-file with the protein pore surface for the last frame.\n")
    
    with hole2.HoleAnalysis(u_last_frame, select='protein', cvect=channel_vec, executable='~/hole2/exe/hole') as ha_last_frame:
        ha_last_frame.run()
        ha_last_frame.create_vmd_surface(filename= new_path+"pore_surface_last_frame.vmd")
        
    ## By using HoleAnalysis as a context manager as shown above,
    ## temporary files are deleted when you are finished with the context manager.



###############################################################################################################################################

## Python main function is the starting point of the program.
## When the program is running, the python interpreter runs the code sequentially.
## In the main function we define the parameters we need for the 'contact_map' function
## by saving the user input into corresponding variables.

def main():

    warnings.filterwarnings("ignore")

    print("\n***** This is a script for analysing pore dimensions of claudin proteins with HOLE2. *****\n")
    print("If the current directory already contains your input files, please press 'Enter' to continue.")
    path_to_dir = input("Otherwise please enter the path to the directory where your input files are saved (e.g. /path/to/input_files):\n   ")
    
    if path_to_dir != "":
        os.chdir(path_to_dir)       ## change directory path

    else:
        path_to_dir = os.getcwd()   ## get the current path

    
    print('\n\nFor protein pore analysis please use input files only containing PROTEIN information.')
    print('These files can be created in VMD by selecting \'protein\' atoms when saving coordinate files during file convertion.')



    top_file = input("\nPlease enter the name of your protein topology file (in DMS format):   ")
    traj_file = input("\nPlease enter the name of your protein trajectory file (in TRR format):   ")
    
    frames = input("\nPlease enter the frame range you want to analyse as a tuple (e.g. (1,101) or (80,101):   ")
    frames = literal_eval(frames)   ## frame range tuple
    
    channel_ax = input("\nAlong which axis is the ion channel in your protein model aligned?\n Please enter one of these options - 'x'/ 'y'/ 'z'):   ")
    
 
    u = mda.Universe(top_file, traj_file, continuous=True)
    u_last_frame = mda.Universe(top_file)
    
    
    ## adding topology attributes to our universe which are not included in DMS-files
        
    u.add_TopologyAttr('record_type')
    u.add_TopologyAttr('altLoc')
    u.add_TopologyAttr('icode')
    u.add_TopologyAttr('occupancy')
    u.add_TopologyAttr('tempfactor')
    u.add_TopologyAttr('element')
    
    u_last_frame.add_TopologyAttr('record_type')
    u_last_frame.add_TopologyAttr('altLoc')
    u_last_frame.add_TopologyAttr('icode')
    u_last_frame.add_TopologyAttr('occupancy')
    u_last_frame.add_TopologyAttr('tempfactor')
    u_last_frame.add_TopologyAttr('element')
    
    
    new_folder_name = input('\nPlease enter a name for a new folder to be created where your output PNG- and CSV-files will be saved:    ')
    os.mkdir(new_folder_name)
    new_path = path_to_dir+'/'+new_folder_name+'/'
    
    print('\n\n************************************************************\n')

    
    hole_analysis(u, u_last_frame, frames, channel_ax, new_path)      ## run the HOLE analysis for the given input variables
    
    
if __name__ == "__main__":
    main()
    
    sys.exit(0)
    
