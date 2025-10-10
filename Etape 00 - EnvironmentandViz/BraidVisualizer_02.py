# Code for thevisualization of the ideal braid. Will start with matplotlib before potentially trying plotly

import numpy as np
import pyvista as pv
import re
from pathlib import Path

# Definitions and parameters for this code.
# Braid details
n = 37 # number of fibers per yarn (from the Braid_Class code)
# Choose plot limits
x_min,x_max = -0.025,0.025
y_min,y_max = -0.025,0.025
z_min,z_max = -0.01,0.06

# Create function for grouping files by yarn type
def extract_letter(file_path):
    match = re.search(r'it_0_(.*?)_fiber',file_path.stem)
    if match:
        return match.group(1)
    else:
        return "" 

# Create function for extracting the number - later they will be put in numeric order
def extract_number(file_path):
    match = re.search(r'(\d+)$',file_path.stem)
    if match:
        return int(match.group(1))
    else:
        return -1

# Prep the plot
plotter = pv.Plotter()

# Define file path and file or files that need to be plotted
archive_dir = Path('c:/Temp/Fiber_Points_step2')
txt_files = sorted(archive_dir.glob('it_0_*_fiber_*.txt'), key=lambda f: (extract_letter(f), extract_number(f)))

# Fiber colours
colours = ['red','blue','green','orange','purple','black']

subset_files = txt_files[:] #+ txt_files[444:444+7*24] # uncomment the second half if there is a second range you want to show

for i,file in enumerate(subset_files):
    with open(file,'r') as f:
        lines = f.readlines()
        # Convert to list of floats then into an array for data processing
        data = [list(map(float,line.strip().split())) for line in lines]
        data = np.array(data)

        # Create mask for filtering the data in the plot so that it is only shown in the desired ranges
        mask = ((data[:,0] >= x_min) & (data[:,0] <= x_max) &
            (data[:,1] >= y_min) & (data[:,1] <= y_max) &
            (data[:,2] >= z_min) & (data[:,2] <= z_max) )
        # and use the created mask
        filtered_data = data[mask]
        
        line = pv.lines_from_points(filtered_data)

        tube = line.tube(radius=0.000150)
        
        # Select the colour based on the colour index for every group of n fibers
        group_index = i // n
        colour = colours[group_index % len(colours)]
        plotter.add_mesh(tube, color=colour)

        #print(f"Plotting file {i+1}: {file.name} with colour '{colour}'") # Useful for debugging the plotting and checking the files

# Position camera as desired
plotter.enable_parallel_projection()  # type: ignore
plotter.camera_position = 'xz'
plotter.window_size = [700,700]

# Control the axis widget and how it looks
plotter.show_axes()  #type: ignore
#plotter.add_axes_at_origin()

# Show the result and potentially take a screenshot, (can also take a screenshot without showing)
desired_dir = 'M:/Recherche/LABSFCA/etudiants/Evans Frandsen/PicturesandSchematics'
plotter.show(
#screenshot = 'desired_dir\Changethis.tiff'
)

#print(plotter.camera_position)