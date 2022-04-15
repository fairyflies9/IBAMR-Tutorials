%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ultimate goal is for this code to write .vertex, .beam, spring, and
% .target files
% Put down dimensions, and calculate dx, ds, and other needed measurements.
Lx = 1;                              % length of computational domain in the x-direction (m)
Ly = 2;                             % length of computational domain in the y-direction (m)
N = 512;                            % number of Cartesian grid meshwidths in the longest direction
dx = Ly/N;                           % Cartesian mesh width (m)
ds = Ly/(2*N);                      % ds should be 1/2 of dx
% Warning, dx = Lx/(N/2) should equal original dx.

% I will have two beams that run horizontally and are Lx/2 in length
Lb = Lx/2;                          %Length of each beam
npts_beam = ceil(Lb/ds);                  %total number of vertices in each beam
npts_total = 2*npts_beam;           %total number of vertices

mesh_name = 'two_beams_';          % structure name
kappa_spring = 2.0e3;               % spring constant (Newton)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make the .vertex file
% Command below is going to create / open a file named mesh_name_N.vertex
vertex_fid = fopen([mesh_name num2str(N) '.vertex'], 'w');

%First line of the file needs to have the total number of vertices
fprintf(vertex_fid, '%d\n', npts_total);

%write out vertices for lower beam
for s = 0:npts_beam-1
   X(1) = Lx/4+s*ds;
   X(2) = Ly/2-Ly/8;
   fprintf(vertex_fid, '%1.16e %1.16e\n', X(1), X(2));
end
disp(['The vertex ids for the bottom beam go from 1 to ' num2str(npts_beam)]);

%write out vertices for upper beam
for s = 0:npts_beam-1
   X(1) = Lx/4+s*ds;
   X(2) = Ly/2+Ly/8;
   fprintf(vertex_fid, '%1.16e %1.16e\n', X(1), X(2));
end
disp(['The vertex ids for the upper beam go from ' num2str(npts_beam+1) ' to ' num2str(npts_total)]);

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%