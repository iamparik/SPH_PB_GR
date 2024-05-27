import re
import math

input_file_path_SA = "bulk.dat"
output_file_path_SA = "input_param_real_particle.dat"


# Regular expressions for GRID and CQUAD4 lines
grid_pattern_SA = re.compile(r'GRID,([\d.-]+),\d+,([Ee\d.+-]+),([Ee\d.+-]+),([Ee\d.+-]+),\d+')
cquad4_pattern = re.compile(r'CQUAD4,([\d.-]+),\d+,([\d.-]+),([\d.-]+),([\d.-]+),([\d.-]+)')

# Extract GRID and CQUAD4 elements
grid_elements_SA = []
cquad4_elements = []

with open(input_file_path_SA, 'r') as input_file:
    for line in input_file:
        grid_match_SA = grid_pattern_SA.match(line)
        cquad4_match = cquad4_pattern.match(line)
        
        if grid_match_SA:
            grid_elements_SA.append(grid_match_SA.groups())
        if cquad4_match:
            cquad4_elements.append(cquad4_match.groups())


# Convert strings to float and store in new lists
meshNodes_SA = []
cquad4Elements= []

for element in grid_elements_SA:
    meshNodes_SA.append(tuple(float(value) for value in element))

for element in cquad4_elements:
    cquad4Elements.append(tuple(float(value) for value in element))

if len(meshNodes_SA) != round(meshNodes_SA[-1][0]):
    print("All nodal points not read")
else:
    print("Total number of nodes in the data file : ", round(meshNodes_SA[-1][0]))

print('size of grid list is', len(meshNodes_SA), "with last element of grid list as",meshNodes_SA)
print('length of cquad is', len(cquad4Elements) , "with cquad list as", cquad4Elements )


# Function to calculate centroid of a QuadElement
def calculate_2Dcentroid(vertices):
    x_SPH = sum(xCAD for _, xCAD, _ , _ in vertices)
    y_SPH = sum(zCAD for _, _, _ , zCAD in vertices)
    num_vertices=float(len(vertices))
    return x_SPH/num_vertices, y_SPH/num_vertices
    


def calculate_area(vertices):
    (_, xA, _, yA), (_, xB, _, yB), (_, xC, _, yC), (_, xD, _, yD)=vertices
    QArea = 0.5 * abs(xA * yB + xB * yC + xC * yD + xD * yA - (yA * xB + yB * xC + yC * xD + yD * xA))
    return QArea

# Calculate centroids of QuadElements and store in a list
centroids = []
QAreas=[]

#print('first element of grid list is', float(grid_elements[0][0]), 'whose data type is', type(float(grid_elements[0][0])))
# Write elements to the output file
with open(output_file_path_SA, 'w') as output_file:

    for quad_element in cquad4Elements:
        Qnumber = 3 #int(quad_element[0])
        vertex_indices = [int(node) for node in quad_element[1:]]
        vertices = [meshNodes_SA[index - 1] for index in vertex_indices]
        centroid = calculate_2Dcentroid(vertices)
        centroids.append(centroid)
        area = calculate_area(vertices)
        QAreas.append(area)
        output_file.write(f"{Qnumber:<10} {centroid[0]/1000.0:<20.10f} {centroid[1]/1000.0:<20.10f} {area/1000000.0:<20.10f}\n")

print("Total Area of All elements : ", sum(QAreas))

print("Extraction and writing of Surface Area Elements completed.")

######### Bdry Calc begin ############################################################################################################

input_file_path_bdry = "bulkbdry.dat"
output_file_path_bdry = "input_param_edge.dat"


# Regular expressions for GRID and CQUAD4 lines
grid_pattern_oneD = re.compile(r'GRID,([\d.-]+),\d+,([Ee\d.+-]+),([Ee\d.+-]+),([Ee\d.+-]+),\d+')
oneD_pattern = re.compile(r'PLOTEL,([\d.-]+),([\d.-]+),([\d.-]+)')

# Extract GRID and CQUAD4 elements
grid_elements_oneD = []
oneD_elements = []

with open(input_file_path_bdry, 'r') as input_file:
    for line in input_file:
        grid_match_oneD = grid_pattern_oneD.match(line)
        oneD_match = oneD_pattern.match(line)
        
        if grid_match_oneD:
            grid_elements_oneD.append(grid_match_oneD.groups())
        if oneD_match:
            oneD_elements.append(oneD_match.groups())


# Convert strings to float and store in new lists
meshNodes_oneD = []
oneDElements= []

for element in grid_elements_oneD:
    meshNodes_oneD.append(tuple(float(value) for value in element))

for element in oneD_elements:
    oneDElements.append(tuple(float(value) for value in element))

if len(meshNodes_oneD) != round(meshNodes_oneD[-1][0]):
    print("All nodal points not read")
else:
    print("Total number of nodes in the data file : ", round(meshNodes_oneD[-1][0]))

print('size of 1D grid list is', len(meshNodes_oneD), "with last element of 1D grid list as",meshNodes_oneD)
print('length of oneDElements is', len(oneDElements) , "with oneD list as", oneDElements )


# Function to calculate centroid of a QuadElement
def calculate_1DmidPoint(vertices):
    x_SPH = sum(xCAD for _, xCAD, _ , _ in vertices)
    y_SPH = sum(zCAD for _, _, _ , zCAD in vertices)
    num_vertices=float(len(vertices))
    return x_SPH/num_vertices, y_SPH/num_vertices
    


def calculate_len_oneD(vertices):
    (_, xA, _, yA), (_, xB, _, yB)=vertices
    distance = math.sqrt((xA- xB)**2 + (yA-yB)**2)
    return distance

# Calculate centroids of QuadElements and store in a list
midPoints = []
oneDlens=[]

#print('first element of grid list is', float(grid_elements[0][0]), 'whose data type is', type(float(grid_elements[0][0])))
# Write elements to the output file
with open(output_file_path_bdry, 'w') as output_file:

    for oneD_element in oneDElements:
        edgeType = 2 #int(quad_element[0])
        vertex_indices = [int(node) for node in oneD_element[1:]]
        vertices = [meshNodes_oneD[index - 1] for index in vertex_indices]
        midPoint = calculate_1DmidPoint(vertices)
        midPoints.append(midPoint)
        oneDlen = calculate_len_oneD(vertices)
        oneDlens.append(oneDlen)
        #Below accounts for XZ plane
        output_file.write(f"{edgeType:<10} {vertices[0][1]/1000.0:<20.10f} {vertices[0][3]/1000.0:<20.10f} {vertices[1][1]/1000.0:<20.10f} {vertices[1][3]/1000.0:<20.10f}\n")

print("Total Area of All elements : ", sum(oneDlens))


######### Create Sim File ############################################################################################################

output_file_path_simSize = "input_param_SimSize.dat"

with open(output_file_path_simSize, 'w') as output_file:
    maxedge=round(len(oneDElements))
    maxnv=round(maxedge+2*len(meshNodes_oneD))
    maxn= maxnv+ len(meshNodes_SA)#len(centroids)
    output_file.write(f"{maxn:<10} {maxnv:<10} {maxedge:<10} {0:<10}\n")
