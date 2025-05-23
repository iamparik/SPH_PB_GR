
#****************************************************************************
#
#  MODULE: NXtoParticle
#
#  PURPOSE:  Uses NX 2D boudary mesh file to generate particles inside the boundary and to 
#           export boundary input_bulkBdryEdge.dat and input_particles_cartesianMesh.dat
#           which can be read by BISPH and BIPI code bases
#
#   CREATED:        08/10/2023       by  PARIKSHIT BOREGOWDA
#   Last Modified:  12/26/2024       by  PARIKSHIT BOREGOWDA 
#****************************************************************************

import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from concurrent.futures import ThreadPoolExecutor


# --- First define all functions to be used ----
# Function to calculate centroid of a TriaElement
def calculate_2Dcentroid(vertices):
    x_SPH = sum(xCAD for _, xCAD, _ , _ in vertices)
    y_SPH = sum(zCAD for _, _, _ , zCAD in vertices)
    num_vertices=float(len(vertices))
    return x_SPH/num_vertices, y_SPH/num_vertices

#Function to area of triangular Element    
def calculate_Tri_area(vertices):
    (_, xA, _, yA), (_, xB, _, yB), (_, xC, _, yC)=vertices
    QArea = 0.5 * abs(xA * (yB - yC) + xB * (yC - yA) + xC * (yA - yB))
    return QArea

#Function to calculate area of Quadrilateral Element
def calculate_Quad_area(vertices):
    (_, xA, _, yA), (_, xB, _, yB), (_, xC, _, yC), (_, xD, _, yD)=vertices
    QArea = 0.5 * abs(xA * yB + xB * yC + xC * yD + xD * yA - (yA * xB + yB * xC + yC * xD + yD * xA))
    return QArea

def calculate_1DmidPoint(vertices):
    x_SPH = sum(xCAD for _, xCAD, _ , _ in vertices)
    y_SPH = sum(zCAD for _, _, _ , zCAD in vertices)
    num_vertices=float(len(vertices))
    return x_SPH/num_vertices, y_SPH/num_vertices
    
#Function to Calculate length of the element
def calculate_len_oneD(vertices):
    (_, xA, _, yA), (_, xB, _, yB)=vertices
    distance = math.sqrt((xA- xB)**2 + (yA-yB)**2)
    return distance

#An algorithm to find if a point is inside polygon for 2D geoemtry
def is_point_inside_polygon(point, polygon,bdry_tolerance):
    x, y = point
    n = len(polygon)
    inside = False

    for i in range(n):
        x1, y1 = polygon[i][0:2]
        x2, y2 = polygon[i][2:4]

        if ((y1 < y and y2 >= y) or (y2 < y and y1 >= y)):
            if x1 + (y - y1) / (y2 - y1) * (x2 - x1) < x:
                inside = not inside # ensure points outsie polygon between parallel lines in same sense are reset to false

    for i in range(n):
        x1, y1 = polygon[i][0:2]
        x2, y2 = polygon[i][2:4]

        dist_xy=abs((x2 - x1)*(y1-y) -(y2 - y1)*(x1-x))/math.sqrt((x2 - x1)*(x2 - x1)+(y2 - y1)*(y2 - y1))
        if (y1 < y and y2 >= y) or (y2 < y and y1 >= y):
            if (dist_xy < bdry_tolerance):
                inside = False

    return inside

#Function to apply Point in Polygon algorithm to all points, using parallel programing
def check_points_parallel(points, polygon, bdry_tolerance,  executor):
    inside_points = []
    outside_points = []

    def check_point(point):
        if is_point_inside_polygon(point, polygon, bdry_tolerance):
            inside_points.append(point)
        else:
            outside_points.append(point)

    futures = [executor.submit(check_point, point) for point in points]
    [future.result() for future in futures]

    return inside_points, outside_points

#Function to apply Point in Polygon algorithm to all points, using serial programing
def check_points_serial(points,polygon, bdry_tolerance):
    inside_points = []
    outside_points = []

    for point in points:
        if is_point_inside_polygon(point, polygon, bdry_tolerance):
            inside_points.append(point)
        else:
            outside_points.append(point)

    return inside_points, outside_points

#Extract values from text lines to tuple containing coordinates of grid points,
# and elements formed by those grid points
def extractValues(line, group_element, text_pattern, elementType):

    group_match_text=text_pattern.match(line)

    # if the text pattern matches the input line, append the group_element
    if group_match_text:
       if(elementType == 0):
           str_group_element = group_match_text.groups()
       else:
           str_group_element = (elementType,)+ group_match_text.groups()

       group_element.append(tuple(float(value) for value in str_group_element))

    return group_element

#Takes NX sim file bulkBdry.dat to obtain bdry elements and particles at cartesian grid coordinates
def bulkBdryCase2D(input_file_directory,output_file_directory, dx, bdry_tolerance,unit_conversion, parallelBool):

    input_file_path=input_file_directory+"bulkBdry1D.dat"
    output_bdryfile_path= output_file_directory+"input_bulkBdryEdge.dat"
    output_cartesianfile_path = output_file_directory+"input_particles_cartesianMesh.dat"
    output_image_path= output_file_directory + "bulkBdry.png"

    if os.path.exists(input_file_path):
        print(" ****************************\n bulkBdry.dat file is processed to : \n generate particles on background grid and \n to process boundary elements \n ****************************")
        
        # Regular expressions for GRID and PLOTEL lines
        grid_pattern_oneD_Bulkbdry = re.compile(r'GRID,([\d.-]+),\d+,([Ee\d.+-]+),([Ee\d.+-]+),([Ee\d.+-]+),\d+')
        oneD_pattern_SolidWall_Bulkbdry = re.compile(r'PLOTEL,([\d.-]+),([\d.-]+),([\d.-]+)')


        # Extract GRID and Plotel line elements
        meshNodes_oneD_Bulkbdry = []
        oneDElements_Bulkbdry= []

        with open(input_file_path, 'r') as input_file:
            for line in input_file:

                #Obtain node (coordinate) values of the grid
                meshNodes_oneD_Bulkbdry= extractValues(line,meshNodes_oneD_Bulkbdry,grid_pattern_oneD_Bulkbdry, 0)
            
                #Find a line element, and identify the nodes of the element
                oneDElements_Bulkbdry= extractValues(line,oneDElements_Bulkbdry,oneD_pattern_SolidWall_Bulkbdry, 2)


        if len(meshNodes_oneD_Bulkbdry) != round(meshNodes_oneD_Bulkbdry[-1][0]):
            print("All 1D nodal points not read")
        else:
            print("Total number of nodes in the data file : ", round(meshNodes_oneD_Bulkbdry[-1][0]))


        # Calculate centroids of QuadElements and store in a list
        midPoints_Bulkbdry = []
        oneDlens_Bulkbdry=[]

        # initialize and define the xmin and y min values
        x_min=[meshNodes_oneD_Bulkbdry[0][1], meshNodes_oneD_Bulkbdry[0][3]]
        x_max=[meshNodes_oneD_Bulkbdry[0][1], meshNodes_oneD_Bulkbdry[0][3]]

        #define bdry polygon
        bdryPolygon=[]

        plt.figure()

        #print('first element of grid list is', float(grid_elements[0][0]), 'whose data type is', type(float(grid_elements[0][0])))
        # Write elements to the output file
        with open(output_bdryfile_path, 'w') as output_file:

            for oneD_element in oneDElements_Bulkbdry:
                vertex_indices = [int(node) for node in oneD_element[2:]]
                vertices = [meshNodes_oneD_Bulkbdry[index - 1] for index in vertex_indices]
                midPoint = calculate_1DmidPoint(vertices)
                midPoints_Bulkbdry.append(midPoint)
                oneDlen = calculate_len_oneD(vertices)
                oneDlens_Bulkbdry.append(oneDlen)
                # Calculate the minimum values in each component
                x_min[0]= min(x_min[0],min(vertices[0][1],vertices[1][1]))
                x_min[1]= min(x_min[1],min(vertices[0][3],vertices[1][3]))
                # Calculate the maximum values in each component
                x_max[0]= max(x_max[0],max(vertices[0][1],vertices[1][1]))
                x_max[1]= max(x_max[1],max(vertices[0][3],vertices[1][3]))
                #create the bdry polygon in he XZ plane
                bdryPolygon.append((vertices[0][1],vertices[0][3],vertices[1][1],vertices[1][3]))
                plt.plot([vertices[0][1],vertices[1][1]],[vertices[0][3],vertices[1][3]], color=(1.0, 0.0, 0.0))
                output_file.write(f"{int(oneD_element[0]):<10} {vertices[0][1]/unit_conversion:<20.10f} {vertices[0][3]/unit_conversion:<20.10f} {vertices[1][1]/unit_conversion:<20.10f} {vertices[1][3]/unit_conversion:<20.10f}\n")

        print("Total Length of All elements : ", sum(oneDlens_Bulkbdry))


        plt.text(0.2, 0.9, f'dx_r = {dx}', fontsize=10, transform=plt.gcf().transFigure)
        plt.text(0.8, 0.9, f'unit {int(unit_conversion)} : 1', fontsize=10, transform=plt.gcf().transFigure)

        plt.savefig(output_image_path) 

        total1DLength_Bulkbdry=float(sum(oneDlens_Bulkbdry))
        totalNum1DElements_Bulkbdry=len(oneDElements_Bulkbdry)
        print("Total length Occupied by _Bulkbdry",totalNum1DElements_Bulkbdry," (all) elements : ", total1DLength_Bulkbdry)
        print("Extraction and writing of Surface Area for _Bulkbdry Elements completed.")
        #------------------------------------------------------------------------------------------

        #------------------------------------------------------------------------------------------
        ############## Generate cartesian cordinates of particles ##########################################   
        totalNum2DElements_Cartesian = fillDomainWParticle2D(x_min, x_max, dx, bdryPolygon, output_cartesianfile_path, bdry_tolerance,unit_conversion, parallelBool)

    else:
        total1DLength_Bulkbdry=0
        totalNum1DElements_Bulkbdry = 0
        totalNum2DElements_Cartesian = 0
    
    return total1DLength_Bulkbdry, totalNum2DElements_Cartesian

# Create cartesian particles for a given boudary information
def fillDomainWParticle2D(x_min, x_max, dx,bdryPolygon, output_file_path, bdry_tolerance,unit_conversion,  parallelBool):
    # Determine size of domain
    size_cartesian = [abs(x_max[0] - x_min[0]),abs(x_max[1] - x_min[1])]

    # Determine grid size
    n_x = math.ceil(size_cartesian[0] / dx) + 1
    n_y = math.ceil(size_cartesian[1] / dx) + 1

    # Find Cartesian points
    cartesianPoints = [[x_min[0]+float(i) * dx +dx/2.0, x_min[1]+float(j) * dx+ dx/2.0] for i in range(0, n_x ) for j in range(0, n_y )]

    # Separate points inside and outside the polygon 
    if(parallelBool):
        # Separate points inside and outside the polygon (in parallel)
        with ThreadPoolExecutor() as executor:
            inside_points, outside_points = check_points_parallel(cartesianPoints, bdryPolygon,bdry_tolerance, executor)
    else:
        # Separate points inside and outside the polygon (in serial)
        inside_points, outside_points = check_points_serial(cartesianPoints, bdryPolygon, bdry_tolerance)

    #Calculate total points
    totalNum2DElements_Cartesian=len(inside_points)

    #Calculate averageArea of elements for the cartesian system
    avgAreaCartesian= dx**2 #total2DArea/totalNum2DElements_Cartesian

    with open(output_file_path, 'w') as output_file:
        for points in inside_points:
            Qnumber = 3 #This is the paritcle type
            output_file.write(f"{Qnumber:<10} {points[0]/unit_conversion:<20.10f} {points[1]/unit_conversion:<20.10f} {avgAreaCartesian/(unit_conversion**2):<20.10f}\n")


    return  totalNum2DElements_Cartesian

#Define the main function
def main():
    #root_directory="../../../"
    #ref_directory="examples/Case9/"  #"examples/Case2/" #"src/Python/NXtoParticle/" 
    #print("current ref_directory is in: ", ref_directory)        

    input_directory="inputCAD/"
    output_directory="data_geo_config/"

    #print("The current code will run from the directory : \n",ref_directory)
    

    # Approximate value of particle spacing
    dx_r = float(input("Enter particle spacing (size of grid) in the .dat file units \n"))

    #Converting Units to meters
    unit_conversion = float(input("\nTo convert input units to new output units, \n enter for 1 output units = ??? input units \n"))

    #tolerance to cutoff a grid point wrt bdry edge
    bdry_tol=dx_r/8.0

    #Confirm parameters used:
    print("\ndx_r =", dx_r, ", \n1 output units = ",unit_conversion, "input units \n")

    #------------------------------------------------------------------------------------------
    ############## Generate 1D edges for the bulk media, using the boundary of the bulk input ##########################################

    input_file_path_Bulkbdry = input_directory
    output_file_path_Bulkbdry = output_directory
    prlBool = True
    totalLen1DElements_Bulkbdry, totalNum2DElements_Cartesian = bulkBdryCase2D(input_file_path_Bulkbdry,output_file_path_Bulkbdry, dx_r, bdry_tol,unit_conversion, prlBool)

    #------------------------------------------------------------------------------------------
    ######### Create Sim File ############################################################################################################

    output_file_path_simSize = output_directory+"/input_param_SimSize.dat"


    with open(output_file_path_simSize, 'w') as output_file:
        output_file.write(f"0 {totalNum2DElements_Cartesian:<10} {totalLen1DElements_Bulkbdry/unit_conversion:<20.30f} 0 0 \n")

    a=input("---------------press return to exit-----------------")

#Now call the main function
if __name__ == "__main__":
    main()
