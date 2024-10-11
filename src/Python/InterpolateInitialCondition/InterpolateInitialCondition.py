import numpy as np
from scipy.interpolate import griddata


def load_semi_analytical_data(file_path):
    """
    Load the data from the given file path.
    The file should have three columns: x, y, and f(x, y).
    """
    data = np.loadtxt(file_path,delimiter =',',dtype=float)
    x = data[:, 0]  # First column: x coordinates
    y = data[:, 1]  # Second column: y coordinates
    f = data[:, 2]  # Third column: function values f(x, y)
    return x, y, f

def load_SPH_coordinate_data(file_path):
    """
    Load the data from the given file path.
    The file should have three columns: x, y, and f(x, y).
    """
    data = np.loadtxt(file_path,dtype=float)
    a = data[:, 1]  # Second column: x coordinates
    b = data[:, 2]  # Third column: y coordinates
    return a,b


def interpolate_value(x, y, f, a, b):
    """
    Interpolate the function value at coordinates (a, b).
    """
    # Create a list of (x, y) points from the data
    points = np.column_stack((x, y))  # Combine x and y into an array of (x, y) pairs

    # Points for interpolation are combined into (a, b) pairs
    grid_points = np.column_stack((a, b))

    # Use griddata to interpolate f at the given (a, b) coordinates
    g = griddata(points, f, grid_points, method='linear')  # You can also use 'nearest' or 'cubic' methods
    return g




#Define the main function
def main():
    root_directory="../../../"
    ref_directory="examples/Case3/"  #"examples/Case3/" #"src/Python/NXtoParticle/" 
    print("current ref_directory is in: ", ref_directory)        
    input_file_path_semiAnalytical = root_directory+ref_directory+"Matlab_PDE_solver/hydroStaticPressure.dat"
    input_file_path_SPH =  root_directory+ref_directory+ "data_geo_config/"+"input_particles_cartesianMesh.dat" # change input_particles_cartesianMesh.dat input_PP.dat to used packed particle positions
    output_file_path=root_directory+ref_directory+ "data_geo_config/"+ "input_particles_param_value.dat"

    print("The current code will run from the directory : \n",ref_directory)

    x, y, f = load_semi_analytical_data(input_file_path_semiAnalytical)
    a,b = load_SPH_coordinate_data(input_file_path_SPH)

    g = interpolate_value(x, y, f, a, b)

    # Write (a, b, g) to a .dat file
    output_data = np.column_stack((a, b, g))
    np.savetxt(output_file_path, output_data,delimiter='\t', fmt='%.6f')
    

#Now call the main function
if __name__ == "__main__":
    main()