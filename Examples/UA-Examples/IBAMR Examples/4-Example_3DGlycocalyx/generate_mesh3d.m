
%This file will generate a rectangular "leaf" that resists bending and stretching. 
%It will be held in place by a line of target points along its center

L = 1;                              % height of computational domain (m)
N = 512;                            % number of Cartesian grid meshwidths at the finest level of the AMR grid
dx = L/N;                           % Cartesian mesh width (m)

glyco_length = 10/64;                  % glyco length (m)
glyco_radius = 1/64;                % glyco radius (m)
glyco_radius = glyco_radius - 1.25*dx;  %glyco radius adjustment for width of delta functions
npts_height = ceil(2*(glyco_length/L)*N);  % number of points along the height
npts_circum = ceil(2*(2*pi*glyco_radius/L)*N); %number of points along the circumference
npts = npts_height*npts_circum;	    % total number points each leaf 
dq = glyco_length/(npts_height-1);     % glyco mesh spacing (m)
ds = glyco_radius/(npts_circum-1);   % glyco mesh spacing (m)

mesh_name = 'glyco_';      % structure name
offset = -0.5*L+dq;	                    % glyco base from bottom of domain (m)

kappa_spring = 5e0;               % spring constant (Newton)
kappa_beam = 0.0; %5.0e-6;                % beam stiffness constant (Newton m^2)
kappa_target = kappa_spring;        % target point penalty spring constant (Newton)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the vertex information

vertex_fid = fopen([mesh_name num2str(N) '.vertex'], 'w');

fprintf(vertex_fid, '%d\n', npts);

for q = 0:npts_height-1
for s = 0:npts_circum-1
   X(1) = 0.0 + glyco_radius*cos(s*2*pi/npts_circum);
   X(2) = 0.0 + glyco_radius*sin(s*2*pi/npts_circum);
   X(3) = offset + q*dq;
   fprintf(vertex_fid, '%1.16e %1.16e %1.16e\n', X(1), X(2), X(3));
end
end 

fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Write out the spring information
% spring_fid = fopen([mesh_name num2str(N) '.spring'], 'w');
% 
% fprintf(spring_fid, '%d\n', npts-npts_height+npts-npts_circum);
% 
% for q = 0:npts_height-1
% for s = 0:npts_circum-2
%    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', q*npts_circum+s, q*npts_circum+s+1, kappa_spring*ds/(ds^2), ds);
% end
% end
% 
% 
% for q = 0:npts_height-2
% for s = 0:npts_circum-1
%    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', q*npts_circum+s, (q+1)*npts_circum+s, kappa_spring*ds/(ds^2), ds);
% end
% end
% 
% fclose(spring_fid);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Write out the beam information
% beam_fid = fopen([mesh_name num2str(N) '.beam'], 'w');
% 
% fprintf(beam_fid, '%d\n', npts-2*npts_height+npts-2*npts_circum);
% 
% for q = 0:npts_height-1
% for s = 0:npts_circum-3
%    fprintf(beam_fid, '%d %d %d %1.16e\n', q*npts_circum+s, q*npts_circum+s+1, q*npts_circum+s+2, kappa_beam*ds/(ds^4));
% end
% end
% 
% 
% for q = 0:npts_height-3
% for s = 0:npts_circum-1
%    fprintf(beam_fid, '%d %d %d %1.16e\n', q*npts_circum+s, (q+1)*npts_circum+s, (q+2)*npts_circum+s, kappa_beam*ds/(ds^4));
% end
% end
% 
% fclose(beam_fid);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Write out the target point information
target_fid = fopen([mesh_name num2str(N) '.target'], 'w');

fprintf(target_fid, '%d\n', npts_height*npts_circum);  %note the number is only npts_chord since we only want a line of target points

for q=0:npts_height-1 %target points will only be attached along a central line
  for s = 0:npts_circum-1
   fprintf(target_fid, '%d %1.16e\n', q*npts_circum+s, kappa_target*ds/(ds^2));
end
end

fclose(target_fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
