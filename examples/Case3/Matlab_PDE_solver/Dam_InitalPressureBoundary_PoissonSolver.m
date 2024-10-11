clear;
global rect_height;
rect_height=input("rectangle height : ");
rect_width=input("rectangle width : ");
Triangle_mesh_Size=input("Triangle_mesh_Size : ");

function bc = HydroP_bc(location,state)
    global rect_height;
    bc = 1000*9.81*(rect_height-location.y);
end

% for rectangle, 3 and 4 are first 2 rows, followed by  4 'x'
% coordiantes and then corresponding 4 'y' cordiantes
 R1=[3,4,0,rect_width,rect_width,0,0,0,rect_height,rect_height]';

%create pde model with 1 pde
model = createpde;

d=decsg(R1);
pdegplot(d,"EdgeLabels","on","FaceLabels","on");



% combine geoemtry and pde with the edges defined
geometryFromEdges(model,d);

% specify coeffecients of the pde - which addresses equations of the form m(u_tt) + d(u_t) −∇⋅(c∇u)+au=f.
specifyCoefficients(model,"m",0,...
                          "d",0,...
                          "c",1,...
                          "a",0,...
                          "f",0);


% Apply the following Neumann boundary conditions on the edge 1.
%  given by : n·(c∇u) + qu = g
applyBoundaryCondition(model,"neumann", ...
                             "Edge",1, ...
                             "g",0, ...
                             "q",0);

%Apply zero Dirichlet condition on the edge 2 and 3.
applyBoundaryCondition(model,"dirichlet", ...
                             "Edge",2:3,"u",0);

% Apply Dirichlet condition h*u = r, where h = 1 and r = rho*g*(y-rect_height)
applyBoundaryCondition(model,"dirichlet", ...
                             "Edge",4, ...
                             "r",@HydroP_bc,"h",1);


% generate mesh
generateMesh(model,"Hmax",Triangle_mesh_Size);
figure
pdemesh(model)

%solve pde
results = solvepde(model);
u = results.NodalSolution;
x= results.Mesh.Nodes(1,:)';
y= results.Mesh.Nodes(2,:)';

figure
pdeplot(model,XYData=u(:,1))

grid_size = length(x);  % Number of grid points in each dimension
xq = linspace(min(x), max(x), grid_size);
yq = linspace(min(y), max(y), grid_size);
[Xq, Yq] = meshgrid(xq, yq);

% Interpolate scattered data to create values on the grid
Fq = griddata(x, y, u, Xq, Yq, 'cubic');  % 'cubic' interpolation

% Plot the contour
figure
contour(Xq, Yq, Fq, 20);  % '20' is the number of contour levels
colorbar;  % Show a colorbar for reference
title('Contour Plot from Scattered Data');
xlabel('X');
ylabel('Y');

% Combine the x, y, and f data into a matrix
data = [x, y, u];

% Get the directory of the current MATLAB file
currentFolder = fileparts(mfilename('fullpath'));

% Create the full path to the .dat file in the same directory as the script
outputFilePath = fullfile(currentFolder, 'hydroStaticPressure.dat');

% Export the data to a .dat file
writematrix(data, outputFilePath, 'Delimiter', ',');


%write hydrostatic dam config used
data= [rect_height, rect_width,Triangle_mesh_Size];
outputFilePath = fullfile(currentFolder, 'hydroStaticPressure_config.dat');
writematrix(data, outputFilePath, 'Delimiter', ';');