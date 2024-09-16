import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

## Example scattered data points (x, y) and their corresponding function values f
#x = np.random.rand(100) * 10  # X coordinates of scattered data
#y = np.random.rand(100) * 10  # Y coordinates of scattered data
#f = np.sin(x) * np.cos(y)  # Example function values f(x, y) based on some function

## New scattered points where you want to interpolate the function values (a, b)
#a = np.random.rand(200) * 10  # X coordinates for interpolation
#b = np.random.rand(200) * 10  # Y coordinates for interpolation

## Interpolating values at the new points (a, b) to get g
## Points for the original data are combined into (x, y) pairs
#points = np.column_stack((x, y))  # Shape: (100, 2)

## Points for interpolation are combined into (a, b) pairs
#grid_points = np.column_stack((a, b))  # Shape: (200, 2)

## Use griddata for interpolation
## Options for method: 'linear', 'nearest', 'cubic'
#g = griddata(points, f, grid_points, method='nearest')

## Compute global min and max for colormap range
#f_min = min(f.min(), g.min())  # Minimum value between both datasets
#f_max = max(f.max(), g.max())  # Maximum value between both datasets


## Plotting the original scattered points and interpolated values for visualization
#plt.figure(figsize=(8, 6))

## Plot original scattered data
#plt.scatter(x, y, c=f, cmap='viridis', label='Original Data', s=50, edgecolor='k', vmin=f_min, vmax=f_max)

## Plot interpolated points with the new values
#plt.scatter(a, b, c=g, cmap='viridis', marker='x', label='Interpolated Data', vmin=f_min, vmax=f_max)

## Add colorbars and labels
#plt.colorbar(label='Function Value')
#plt.xlabel('X')
#plt.ylabel('Y')
#plt.title('Interpolation of Scattered Data in 2D')

## Show the legend and plot
#plt.legend()
#plt.grid(True)
#plt.show()


def load_data(file_path):
    """
    Load the data from the given file path.
    The file should have three columns: x, y, and f(x, y).
    """
    data = np.loadtxt(file_path)
    x = data[:, 0]  # First column: x coordinates
    y = data[:, 1]  # Second column: y coordinates
    f = data[:, 2]  # Third column: function values f(x, y)
    return x, y, f



#Define the main function
def main():
    root_directory="../../../"
    ref_directory="examples/Case3/"  #"examples/Case3/" #"src/Python/NXtoParticle/" 
    print("current ref_directory is in: ", ref_directory)        

    input_file_path="Matlab_PDE_solver/hydroStaticPressure.dat"
    output_directory="data_geo_config/"

    print("The current code will run from the directory : \n",ref_directory)

    x, y, f = load_data(input_file_path)

    print("x = ", x)
    print("y = ", y)
    print("f = ", f)

    

#Now call the main function
if __name__ == "__main__":
    main()