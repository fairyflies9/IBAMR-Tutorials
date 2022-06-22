L = 1;                              % length of computational domain (m)
N = 512;                            % number of Cartesian grid meshwidths at the finest level of the AMR grid
dx = L/N;                           % Cartesian mesh width (m)

center_x = 0;                     %x-coordinate of the center of the circle.
center_y = 0;                     %y-coordinate of the center of the circle.

radius = 0.1;                       % radius of band (m)
epsilon = 0.02;                    %deformation in the x-direction
band_length = 2*pi*radius;          % rubber band length (m)
npts = ceil(2*(band_length/L)*N);   % number of points along the rubber band
ds = band_length/(npts);            % physical distance between neighboring Lagrangian mesh points (m)
dtheta = 2*pi/npts;

mesh_name = 'rubber_band_';          % structure name
kappa_spring = 2.0e3;               % spring constant (Newton)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the vertex information
vertex_fid = fopen([mesh_name num2str(N) '.vertex'], 'w');

fprintf(vertex_fid, '%d\n', npts);

for s = 0:npts-1
   X(1) = (radius+epsilon)*sin(s*dtheta)+center_x;
   X(2) = radius*cos(s*dtheta)+center_y;
   fprintf(vertex_fid, '%1.16e %1.16e\n', X(1), X(2));
end

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the spring information
spring_fid = fopen([mesh_name num2str(N) '.spring'], 'w');

fprintf(spring_fid, '%d\n', npts);

for s = 0:npts-2
   fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, kappa_spring*ds/(ds^2), ds);
end


fprintf(spring_fid, '%d %d %1.16e %1.16e\n', npts-1, 0, kappa_spring*ds/(ds^2), ds);

fclose(spring_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




