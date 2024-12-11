
from scipy.interpolate import griddata
import pandas as pd
import numpy as np
import os
import re

def get_max_file_index(folder_path, file_pattern):
    """
    Find the maximum file index in the specified directory based on the file pattern.
    
    Parameters:
    - folder_path: Path to the folder where the files are stored.
    - file_pattern: Pattern of file names, with '{}' placeholder for the file index.
    
    Returns:
    - max_index: Maximum file index found.
    """
    max_index = 0
    file_regex = re.compile(file_pattern.replace("{}", r"(\d+)"))
    
    for file_name in os.listdir(folder_path):
        match = file_regex.match(file_name)
        if match:
            file_index = int(match.group(1))
            max_index = max(max_index, file_index)
    
    return max_index


def load_SPH_data(file_path):
    """
    Load SPH data from a given file path, skipping the header row.
    
    Parameters:
    - file_path: Path to the file to load.
    
    Returns:
    - df: DataFrame with columns x, y, vx, vy, p if file exists, otherwise None.
    """
    if os.path.exists(file_path):
        return pd.read_csv(file_path, delimiter=',', skiprows=1, names=['x', 'y', 'vx', 'vy', 'p'])
    else:
        print(f"File {file_path} does not exist.")
        return None

def load_CFD_data(file_path):
    """
    Load CFD data from the given file path using pandas.
    Expects a CSV file with columns: x, y, p, vx, vy.
    
    Parameters:
    - file_path: Path to the CFD data file.
    
    Returns:
    - df: DataFrame with standardized column names.
    """
    # Read data with headers, assuming the file has column names in the first row
    df = pd.read_csv(file_path)

    # Rename columns to standardized names
    df.rename(columns={
        "Absolute Pressure (Pa)": "p",
        "Velocity[i] (m/s)": "vx",
        "Velocity[j] (m/s)": "vy",
        "X (m)": "x",
        "Y (m)": "y"
    }, inplace=True)
    return df

def interpolate_value(x, y, f, a, b):
    """
    Interpolate the function value at coordinates (a, b).
    
    Parameters:
    - x, y: Arrays of coordinates for known data points.
    - f: Array of function values at the known data points.
    - a, b: Coordinates where we want to interpolate the function values.
    
    Returns:
    - g: Interpolated values at (a, b).
    """
    # Combine x and y into an array of (x, y) points
    points = np.column_stack((x, y))

    # Combine a and b into an array of (a, b) points for interpolation
    grid_points = np.column_stack((a, b))

    # Use griddata to interpolate f at the given (a, b) coordinates
    g = griddata(points, f, grid_points, method='linear')  # Options: 'linear', 'nearest', 'cubic'
    return g

def calculate_errors(df_CFD, df_SPH):
    """
    Calculate the interpolation error between CFD and SPH data.

    Parameters:
    - df_CFD: DataFrame containing CFD data with columns x, y, vx, vy, p.
    - df_SPH: DataFrame containing SPH data with columns x, y, vx, vy, p.

    Returns:
    - err_p: Mean error for pressure.
    - err_vx: Mean error for velocity in x direction.
    - err_vy: Mean error for velocity in y direction.
    """
    # Interpolate CFD values to SPH points
    p_ref = interpolate_value(df_CFD['x'].values, df_CFD['y'].values, df_CFD['p'].values, df_SPH['x'].values, df_SPH['y'].values)
    vx_ref = interpolate_value(df_CFD['x'].values, df_CFD['y'].values, df_CFD['vx'].values, df_SPH['x'].values, df_SPH['y'].values)
    vy_ref = interpolate_value(df_CFD['x'].values, df_CFD['y'].values, df_CFD['vy'].values, df_SPH['x'].values, df_SPH['y'].values)

    # Calculate mean error for each parameter
    err_p = np.sqrt(((df_SPH['p'] - p_ref) ** 2).mean()) 
    err_vx = np.sqrt(((df_SPH['vx'] - vx_ref) ** 2).mean())
    err_vy = np.sqrt(((df_SPH['vy'] - vy_ref) ** 2).mean()) 
    return err_p, err_vx, err_vy

def save_errors(output_path, errs):
    """
    Save the computed averages to a file.
    
    Parameters:
    - output_path: Path to the output file.
    - averages: List of tuples with average data to save.
    """
    errs_df = pd.DataFrame(errs, columns=['file_index', 'err_p', 'err_vx', 'err_vy'])
    errs_df.to_csv(output_path, sep='\t', index=False, float_format='%.6f')
    print(f"Averages saved to {output_path}")


def process_files(folder_path, file_pattern,df_CFD):
    """
    Process all files in the specified folder based on the file pattern and the maximum file index,
    calculating averages for each file and storing them in a list.
    
    Parameters:
    - folder_path: Path to the folder where the files are stored.
    - file_pattern: Pattern of file names to read (e.g., "p{}.dat").
    
    Returns:
    - List of tuples containing file index and average values (avg_x, avg_y, avg_vx, avg_vy, avg_p).
    """
    errors = []
    max_index = get_max_file_index(folder_path, file_pattern)
    
    for i in range(1, max_index + 1):
        file_path = os.path.join(folder_path, file_pattern.format(i))
        df_SPH = load_SPH_data(file_path)
        
        if df_SPH is not None:
            print("currently in file number ", i)
            err_p, err_vx, err_vy = calculate_errors(df_CFD, df_SPH)
            errors.append((i, err_p, err_vx, err_vy))
    
    return errors




#Define the main function
def main():
    folder_path = "../data_output"  # Directory containing the data files
    file_pattern_SPH = "p{}.dat"
    input_file_CFD = "../CFD_Data/inititalValuesCFD.csv"
    output_file_path = "error_output.dat"

    df_CFD=load_CFD_data(input_file_CFD)
    
    # Process files and calculate averages
    errs = process_files(folder_path, file_pattern_SPH, df_CFD)
    
    # Save averages to output file
    save_errors(output_file_path, errs)


#Now call the main function
if __name__ == "__main__":
    main()
