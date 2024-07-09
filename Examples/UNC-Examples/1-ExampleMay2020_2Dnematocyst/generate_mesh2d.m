L = 0.0005;                              % length of computational domain (m)
N = 512;                            % number of Cartesian grid meshwidths at the finest level of the AMR grid
dx = L/N;                           % Cartesian mesh width (m)

radius = 0.00005/2;               % radius of band (m)
epsilon = 0.00;                    % deformation in the x-direction
band_length = 2*pi*radius;          % rubber band length (m)
npts = ceil(2*(band_length/L)*N);   % number of points along the rubber band
ds = band_length/(npts);            % physical distance between neighboring Lagrangian mesh points (m)
dtheta = 2*pi/npts;
prey_distance = L/10;               %distance between beam tip and prey. Note L/10 = beam length
offset = -L/4+(0.5)*L/10+prey_distance+radius;

mesh_name = 'rubber_band';          % structure name

kappa_spring = 8.0e4;               % spring constant (Newton)
kappa_beam = 8.9e-9;                % beam stiffness constant (Newton m^2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the vertex information
vertex_fid = fopen([mesh_name num2str(N) '.vertex'], 'w');

fprintf(vertex_fid, '%d\n', npts);

for s = 0:npts-1
   X(1) = (radius+epsilon)*sin(s*dtheta);
   X(2) = offset+radius*cos(s*dtheta);
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

% Write out the beam information

beam_fid = fopen([mesh_name num2str(N) '.beam'], 'w');

fprintf(beam_fid, '%d\n', npts);

for s = 0:npts-3
   fprintf(beam_fid, '%d %d %d %1.16e\n', s, s+1, s+2, kappa_beam*ds/(ds^4));
end

fprintf(beam_fid, '%d %d %d %1.16e\n', npts-2, npts-1, 0, kappa_beam*ds/(ds^4));
fprintf(beam_fid, '%d %d %d %1.16e\n', npts-1, 0, 1, kappa_beam*ds/(ds^4));
fclose(beam_fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


