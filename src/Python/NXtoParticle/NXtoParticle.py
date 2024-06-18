#****************************************************************************
#
#  MODULE: NXtoParticle
#
#  PURPOSE:  General Configuration of file is set here
#
#   CREATED:        08/10/2023       by  PARIKSHIT BOREGOWDA
#   Last Modified:  04/13/2024       by  PARIKSHIT BOREGOWDA 
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

#Takes NX sim file bulk.dat and converts elements to particles for SPH code
def bulkCase2D(input_file_directory, output_file_directory,unit_conversion):

    input_file_path= input_file_directory+"bulk2D.dat"
    output_file_path = output_file_directory+"input_particles_bulkMesh.dat"

    if os.path.exists(input_file_path):
        print(" ****************************\n bulk2D.dat file is processed to \n generate particles \n **********************")

        # Regular expressions for GRID, CTRIA3 and CQUAD4 lines
        grid_pattern_SA = re.compile(r'GRID,([\d.-]+),\d+,([Ee\d.+-]+),([Ee\d.+-]+),([Ee\d.+-]+),\d+')
        ctria3_pattern = re.compile(r'CTRIA3,([\d.-]+),\d+,([\d.-]+),([\d.-]+),([\d.-]+)')
        cquad4_pattern = re.compile(r'CQUAD4,([\d.-]+),\d+,([\d.-]+),([\d.-]+),([\d.-]+),([\d.-]+)')


        # Extract GRID nodes, TRIA3 and CQUAD4 elements
        grid_nodes = []
        ctria3Elements= []
        cquad4Elements= []

        with open(input_file_path, 'r') as input_file:
            for line in input_file:
                #Obtain node (coordinate) values of the grid
                grid_nodes= extractValues(line,grid_nodes,grid_pattern_SA, 0)

                #Find a triangle element, and identify the nodes of the element
                ctria3Elements= extractValues(line,ctria3Elements,ctria3_pattern, 3)

                #Find a Quadrilateral element, and identify the nodes of the element
                cquad4Elements= extractValues(line,cquad4Elements,cquad4_pattern, 3)

        if len(grid_nodes) != round(grid_nodes[-1][0]):
            print("All area element's nodal points not read")
        else:
            print("Total number of nodes in the data file : ", round(grid_nodes[-1][0]))


        # Calculate centroids of QuadElements and store in a list
        centroids = []
        QAreas=[]

        #print('first element of grid list is', float(grid_elements[0][0]), 'whose data type is', type(float(grid_elements[0][0])))
        # Write elements to the output file
        with open(output_file_path, 'w') as output_file:
    
            for tria_element in ctria3Elements:
                vertex_indices = [int(node) for node in tria_element[2:]]
                vertices = [grid_nodes[index - 1] for index in vertex_indices]
                centroid = calculate_2Dcentroid(vertices)
                centroids.append(centroid)
                area = calculate_Tri_area(vertices)
                QAreas.append(area)
                output_file.write(f"{int(tria_element[0]):<10} {centroid[0]/unit_conversion:<20.10f} {centroid[1]/unit_conversion:<20.10f} {area/(unit_conversion**2):<20.10f}\n")

    
            for quad_element in cquad4Elements:
                vertex_indices = [int(node) for node in quad_element[2:]]
                vertices = [grid_nodes[index - 1] for index in vertex_indices]
                centroid = calculate_2Dcentroid(vertices)
                centroids.append(centroid)
                area = calculate_Quad_area(vertices)
                QAreas.append(area)
                output_file.write(f"{int(quad_element[0]):<10} {centroid[0]/unit_conversion:<20.10f} {centroid[1]/unit_conversion:<20.10f} {area/(unit_conversion**2):<20.10f}\n")

            

        total2DArea=sum(QAreas)
        totalNum2DElements=len(ctria3Elements)+len(cquad4Elements)
        AvgArea=total2DArea/totalNum2DElements
        DiaofParticle=math.sqrt(AvgArea*4.0/math.pi)

        print("Total Area Occupied by",totalNum2DElements," (all) elements : ", total2DArea)
        print("Average Element Area : ", AvgArea ," gives average particle diameter : ", DiaofParticle)

        print("Extraction and writing of Surface Area Elements completed.")
    else :
        totalNum2DElements = 0
        total2DArea = 0


    

    return totalNum2DElements, total2DArea

#Takes NX sim file bulkBdry.dat to obtain bdry elements and particles at cartesian grid coordinates
def bulkBdryCase2D(input_file_directory,output_file_directory, dx, bdry_tolerance,unit_conversion, parallelBool):

    input_file_path=input_file_directory+"bulkBdry1D.dat"
    output_bdryfile_path= output_file_directory+"input_bulkBdryEdge.dat"
    output_cartesianfile_path = output_file_directory+"input_particles_cartesianMesh.dat"

    if os.path.exists(input_file_path):
        print(" ****************************\n bulkBdry.dat file is processed to : \n generate particles on background grid and \n to process boundary elements \n ****************************")
        
        # Regular expressions for GRID and CQUAD4 lines
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
        verticesEdge_Bulkbdry=[]

        # initialize and define the xmin and y min values
        x_min=[meshNodes_oneD_Bulkbdry[0][1], meshNodes_oneD_Bulkbdry[0][3]]
        x_max=[meshNodes_oneD_Bulkbdry[0][1], meshNodes_oneD_Bulkbdry[0][3]]

        #define bdry polygon
        bdryPolygon=[]

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
                output_file.write(f"{int(oneD_element[0]):<10} {vertices[0][1]/unit_conversion:<20.10f} {vertices[0][3]/unit_conversion:<20.10f} {vertices[1][1]/unit_conversion:<20.10f} {vertices[1][3]/unit_conversion:<20.10f}\n")

        print("Total Length of All elements : ", sum(oneDlens_Bulkbdry))

        total1DLength_Bulkbdry=sum(oneDlens_Bulkbdry)
        totalNum1DElements_Bulkbdry=len(oneDElements_Bulkbdry)
        print("Total length Occupied by _Bulkbdry",totalNum1DElements_Bulkbdry," (all) elements : ", total1DLength_Bulkbdry)
        print("Extraction and writing of Surface Area for _Bulkbdry Elements completed.")
        #------------------------------------------------------------------------------------------

        #------------------------------------------------------------------------------------------
        ############## Generate cartesian cordinates of particles ##########################################   
        totalNum2DElements_Cartesian = fillDomainWParticle2D(x_min, x_max, dx, bdryPolygon, output_cartesianfile_path, bdry_tolerance,unit_conversion, parallelBool)

    else:
        totalNum1DElements_Bulkbdry = 0
        totalNum2DElements_Cartesian = 0
    
    return totalNum1DElements_Bulkbdry, totalNum2DElements_Cartesian

#Takes NX sim file domainBdry.dat to obtain bdry elements and particles at cartesian grid coordinates
def domainBdryCase2D(input_file_directory, output_file_directory, dx, unit_conversion):

    input_file_path= input_file_directory +"domainBdry1D.dat"
    output_file_path= output_file_directory +"input_domainBdryEdge.dat"
    output_image_path= output_file_directory + "domainBdry.png"

    if os.path.exists(input_file_path):
        print(" **************************** \n domainBdry.dat file is processed to : \n to process boundary elements of domain \n ****************************")
        # Regular expressions for GRID and CQUAD4 lines
        grid_pattern_oneD = re.compile(r'GRID,([\d.-]+),\d+,([Ee\d.+-]+),([Ee\d.+-]+),([Ee\d.+-]+),\d+')
        oneD_pattern_SolidWall = re.compile(r'PLOTEL,([\d.-]+),([\d.-]+),([\d.-]+)')
   
        # Extract GRID and CTRIA3 elements
        meshNodes_oneD = []
        oneDElements = []

        bdryType = 0
        with open(input_file_path, 'r') as input_file:
            for line in input_file:

                if(re.search(r'Mesh Collector',line)):
                    bdryType=bdryType+1

                #Obtain node (coordinate) values of the grid
                meshNodes_oneD= extractValues(line,meshNodes_oneD,grid_pattern_oneD, 0)

                #Find a line element, and identify the nodes of the element
                oneDElements= extractValues(line,oneDElements,oneD_pattern_SolidWall, bdryType)


        if len(meshNodes_oneD) != round(meshNodes_oneD[-1][0]):
            print("All 1D nodal points not read")
        else:
            print("Total number of nodes in the data file : ", round(meshNodes_oneD[-1][0]))

        #print('size of 1D grid list is', len(meshNodes_oneD), "with last element of 1D grid list as",meshNodes_oneD)
        #print('length of oneDElements is', len(oneDElements) , "with oneD list as", oneDElements )


        # Calculate centroids of QuadElements and store in a list
        midPoints = []
        oneDlens=[]


        plt.figure()
        # Create a dictionary to map integers to colors
        int_to_color = {
            1: (1.0, 0.0, 0.0),    # Red
            2: (0.0, 1.0, 0.0),    # Green
            3: (0.0, 0.0, 1.0),    # Blue
            4: (1.0, 1.0, 0.0),    # Yellow
            5: (1.0, 0.0, 1.0),    # Magenta
            6: (0.0, 1.0, 1.0),    # Cyan
            7: (0.5, 0.5, 0.5),    # Gray
            8: (1.0, 0.5, 0.0),    # Orange
            9: (0.5, 0.0, 0.5),    # Purple
            10: (0.0, 0.5, 0.5),   # Teal
            # Add more mappings as needed
            }
        #print('first element of grid list is', float(grid_elements[0][0]), 'whose data type is', type(float(grid_elements[0][0])))
        # Write elements to the output file
        with open(output_file_path, 'w') as output_file:

            for oneD_element in oneDElements:         
                vertex_indices = [int(node) for node in oneD_element[2:]]
                vertices = [meshNodes_oneD[index - 1] for index in vertex_indices]
                midPoint = calculate_1DmidPoint(vertices)
                midPoints.append(midPoint)
                oneDlen = calculate_len_oneD(vertices)
                oneDlens.append(oneDlen)
                #Below accounts for XZ plane
                output_file.write(f"{int(oneD_element[0]):<10} {vertices[0][1]/unit_conversion:<20.10f} {vertices[0][3]/unit_conversion:<20.10f} {vertices[1][1]/unit_conversion:<20.10f} {vertices[1][3]/unit_conversion:<20.10f}\n")
                color=int_to_color.get(int(oneD_element[0]),(0.0,0.0,0.0)) # Default to black if color is not found
                plt.plot([vertices[0][1],vertices[1][1]],[vertices[0][3],vertices[1][3]], color=color)

        # now we print bdry color labels corresponsing to the assigned unique bdry number
        for label, color in int_to_color.items():
            plt.text(0.2, 0.9, f'dx_r = {dx}', fontsize=10, transform=plt.gcf().transFigure)
            plt.text(0.8, 0.9, f'unit {int(unit_conversion)} : 1', fontsize=10, transform=plt.gcf().transFigure)
            if (label <= bdryType):
                color_name = f'BC {label}'
                plt.text(0.02, 0.95 - label * 0.05, f'{color_name}', fontsize=10,
                     transform=plt.gcf().transFigure, va='top', bbox=dict(facecolor=color, alpha=0.5))


        plt.savefig(output_image_path) 
        #plt.show()

        print("Total Length of All elements : ", sum(oneDlens))

        total1DLength=sum(oneDlens)
        totalNum1DElements=len(oneDElements)
        print("Total length Occupied by",totalNum1DElements," (all) elements : ", total1DLength)
        print("Extraction and writing of Surface Area Elements completed.")
        #------------------------------------------------------------------------------------------
    else:
        totalNum1DElements = 0


    

    return totalNum1DElements

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
    root_directory="../../../"
    ref_directory="src/Python/NXtoParticle/"  #"examples/Case2/" #"src/Python/NXtoParticle/" 
    print("current ref_directory is in: ", ref_directory)        

    input_directory="inputCAD/"
    output_directory="data_geo_config/"

    print("The current code will run from the directory : \n",ref_directory)
    

    # Approximate value of particle spacing
    dx_r = float(input("Enter particle spacing (size of grid) in the .dat file units \n"))

    #Converting Units to meters
    unit_conversion = float(input("\nTo convert input units to new output units, \n enter for 1 output units = ??? input units \n"))

    #tolerance to cutoff a grid point wrt bdry edge
    bdry_tol=dx_r/8.0

    #Confirm parameters used:
    print("\ndx_r =", dx_r, ", \n1 output units = ",unit_conversion, "input units \n")

    
    #------------------------------------------------------------------------------------------
    ############# First convert area mesh generated in CAD to 2D particles #####################
    input_file_path_SA = root_directory+ref_directory+input_directory
    output_file_path_SA = root_directory+ref_directory+output_directory

    totalNum2DElements,total2DArea= bulkCase2D(input_file_path_SA,output_file_path_SA,unit_conversion)
    #------------------------------------------------------------------------------------------

    #------------------------------------------------------------------------------------------
    ############## Generate 1D edges for the bulk media, using the boundary of the bulk input ##########################################

    input_file_path_Bulkbdry = root_directory+ref_directory+input_directory
    output_file_path_Bulkbdry = root_directory+ref_directory+output_directory
    prlBool = True
    totalNum1DElements_Bulkbdry, totalNum2DElements_Cartesian = bulkBdryCase2D(input_file_path_Bulkbdry,output_file_path_Bulkbdry, dx_r, bdry_tol,unit_conversion, prlBool)

        #------------------------------------------------------------------------------------------

    #------------------------------------------------------------------------------------------
    ############## Generate 1D edges for the domain boundary ##########################################

    input_file_path_bdry = root_directory+ref_directory+input_directory
    output_file_path_bdry = root_directory+ref_directory+output_directory

    totalNum1DElements=domainBdryCase2D(input_file_path_bdry,output_file_path_bdry,dx_r, unit_conversion)

    #------------------------------------------------------------------------------------------
    ######### Create Sim File ############################################################################################################

    output_file_path_simSize = root_directory+ref_directory+output_directory+"/input_param_SimSize.dat"


    with open(output_file_path_simSize, 'w') as output_file:
        output_file.write(f"{totalNum2DElements:<10} {totalNum2DElements_Cartesian:<10} {totalNum1DElements_Bulkbdry:<10} {totalNum1DElements:<10} {total2DArea/(unit_conversion**2):<20.30f}\n")

    a=input("---------------press return to exit-----------------")

#Now call the main function
if __name__ == "__main__":
    main()
