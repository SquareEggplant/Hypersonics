function [T, P, U, V, W, centroids, norms, areas] = load_vtk_surf(vtkFile)
% Loads an ASCII VTK POLYDATA surface file into a MATLAB table.
% Each row is a triangle with geometry and all cell data fields.
%
% Output table columns:
%   - area: scalar (triangle area)
%   - normal: 1x3 (unit normal)
%   - centroid: 1x3 (centroid coordinates)
fid = fopen(vtkFile,'r');
assert(fid>0, 'Cannot open file.');
% --- Skip header lines until POINTS section
while true
   line = fgetl(fid);
   if contains(line, 'POINTS')
       break;
   end
end
% --- Read points
tokens = split(line);
points = str2double(tokens(2));
pointarr = zeros(points,3);
for i = 1:points
   line = fgetl(fid);
   pointarr(i,:) = str2double(split(line));
end
% --- Skip header lines until POLYGONS section
while true
   line = fgetl(fid);
   if contains(line, 'POLYGONS')
       break;
   end
end
% --- Read POLYGONS
tokens = split(line);
num_tri = str2double(tokens(2));
%Table declarations
norms = zeros(num_tri, 3);
areas = zeros(num_tri, 1);
P = zeros(num_tri, 1);
T = zeros(num_tri, 1);
U = zeros(num_tri, 1);
V = zeros(num_tri, 1);
W = zeros(num_tri, 1);
centroids = zeros(num_tri, 3);
for i = 1:num_tri
   line = fgetl(fid);
   tokens = str2double(split(line));
   v1 = pointarr(tokens(2)+1,:);
   v2 = pointarr(tokens(3)+1,:);
   v3 = pointarr(tokens(4)+1,:);
   trinorm = cross(v2-v1,v3-v2);
   areas(i) = 0.5*norm(trinorm);
   centroids(i,:) = (v1+v2+v3)/3;
   norms(i,:) = trinorm/norm(trinorm);
   %if norms(i,:)==0
   %    norms(i,:)=[-1,0,0];
   %end
end
%%
% --- Skip header lines until Pressure section
while true
   line = fgetl(fid);
   if contains(line, ' P ')
       break;
   end
end
   line = fgetl(fid);
% --- Read Pressures
for i = 1:num_tri
   line = fgetl(fid);
   P(i) = str2double(line);
end
%%
% --- Skip header lines until Temp section
while true
   line = fgetl(fid);
   if contains(line, ' T ')
       break;
   end
end
   line = fgetl(fid);
% --- Read Temps
for i = 1:num_tri
   line = fgetl(fid);
   T(i) = str2double(line);
end
% --- Skip header lines until U section
while true
   line = fgetl(fid);
   if contains(line, ' U ')
       break;
   end
end
   line = fgetl(fid);
% --- Read U
for i = 1:num_tri
   line = fgetl(fid);
   U(i) = str2double(line);
end
% --- Skip header lines until V section
while true
   line = fgetl(fid);
   if contains(line, ' V ')
       break;
   end
end
   line = fgetl(fid);
% --- Read V
for i = 1:num_tri
   line = fgetl(fid);
   V(i) = str2double(line);
end
% --- Skip header lines until W section
while true
   line = fgetl(fid);
   if contains(line, ' W ')
       break;
   end
end
   line = fgetl(fid);
% --- Read W
for i = 1:num_tri
   line = fgetl(fid);
   W(i) = str2double(line);
end

