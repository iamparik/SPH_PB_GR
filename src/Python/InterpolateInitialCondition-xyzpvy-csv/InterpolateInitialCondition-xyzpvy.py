import numpy as np
from scipy.interpolate import griddata


def load_CFDdata(file_path):
    """
    Load the data from the given file path.
    The file should have three columns: x, y, and f(x, y).
    """
    data = np.loadtxt(file_path, delimiter=',', dtype=float, skiprows=1)
    x = data[:, 5]  # First column: x coordinates
    y = data[:, 6]  # Second column: y coordinates
    p = data[:, 0]  # Third column: pressure values p(x, y)
    vx = data[:, 1]  # Third column: function values vx(x, y)
    vy = data[:, 2]  # Third column: function values vy(x, y)
    return x, y, vx,vy, p

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
    g = griddata(points, f, grid_points, method='linear')  # You can  use 'linear' or 'nearest' or 'cubic' methods
    return g




#Define the main function
def main():
    root_directory="../../../"
    ref_directory="examples/Case1c/"  #"examples/Case3/" #"src/Python/NXtoParticle/" 
    print("current ref_directory is in: ", ref_directory)        
    input_file_path_semiAnalytical = root_directory+ref_directory+"CFD_Data/inititalValuesCFD.csv"
    input_file_path_SPH1 =  root_directory+ref_directory+ "data_geo_config/"+"input_PP.dat" # change input_particles_cartesianMesh.dat as input_PP.dat to use packed particle positions
    input_file_path_SPH2 =  root_directory+ref_directory+ "data_geo_config/"+"input_particles_cartesianMesh.dat"
    output_file_path1=root_directory+ref_directory+ "data_geo_config/"+ "input_PPparticles_param_value.dat"
    output_file_path2=root_directory+ref_directory+ "data_geo_config/"+ "input_Cartparticles_param_value.dat"

    
    print("The current code will run from the directory : \n",ref_directory)

    x, y, vx,vy, p = load_CFDdata(input_file_path_semiAnalytical)
    a1,b1 = load_SPH_coordinate_data(input_file_path_SPH1)
    a2,b2 = load_SPH_coordinate_data(input_file_path_SPH2)

    p_new1 = interpolate_value(x, y, p, a1, b1)
    vx_new1 = interpolate_value(x, y,vx, a1, b1)
    vy_new1 = interpolate_value(x, y,vy, a1, b1)

    p_new2 = interpolate_value(x, y, p, a2, b2)
    vx_new2 = interpolate_value(x, y,vx, a2, b2)
    vy_new2 = interpolate_value(x, y,vy, a2, b2)

    # Write (a, b, g) to a .dat file
    output_data1 = np.column_stack((a1, b1, vx_new1, vy_new1, p_new1))
    output_data2 = np.column_stack((a2, b2, vx_new2, vy_new2, p_new2))
    np.savetxt(output_file_path1, output_data1,delimiter='\t', fmt='%.6f')
    np.savetxt(output_file_path2, output_data2,delimiter='\t', fmt='%.6f')
    

#Now call the main function
if __name__ == "__main__":
    main()