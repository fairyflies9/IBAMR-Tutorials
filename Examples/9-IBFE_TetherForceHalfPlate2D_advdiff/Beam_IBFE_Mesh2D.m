close all;
clear all;
L = 2;                              % length of computational domain (m)
MAX_LEVELS = 2;                            % maximum number of levels in locally refined grid
REF_RATIO  = 4;                           % refinement ratio between levels
NCOARSE = 32;                                     % actual    number of grid cells on coarsest grid level
N = (REF_RATIO^(MAX_LEVELS - 1))*NCOARSE;  % effective number of grid cells on finest   grid level
dx = (1.0*L)/N;                           % Cartesian mesh width (m)
ds = dx*2;


x_points(1:33)=-.5:ds:.5;
x_points(34:66)=-.5:ds:.5;

y_points(1:33)=.125;
y_points(34:66)=-.125;

y_points(67:73)=(-.125+ds):ds:(.125-ds);
x_points(67:73)=.5;

y_points(74:80)=(-.125+ds):ds:(.125-ds);
x_points(74:80)=-.5;

id_points=1:80;

p=80;
%% CONFIGURING EDGES
facets(1:32,1)=id_points(1:32);
facets(1:32,2)=id_points(2:33);

facets(33:64,1)=id_points(34:65);
facets(33:64,2)=id_points(35:66);

facets(65:70,1)=id_points(67:72);
facets(65:70,2)=id_points(68:73);

facets(71:76,1)=id_points(74:79);
facets(71:76,2)=id_points(75:80);

facets(77,1)=id_points(1);
facets(77,2)=id_points(80);
facets(78,1)=id_points(33);
facets(78,2)=id_points(73);
facets(79,1)=id_points(34);
facets(79,2)=id_points(74);
facets(80,1)=id_points(66);
facets(80,2)=id_points(67);


%% CREATING THE INPUT FILES
mesh_name = 'IBFE_Mesh2D_';


vertex_fid = fopen([mesh_name num2str(N) '.poly'], 'w');
fprintf(vertex_fid, '%d %d %d %d \n \n', p, 2, 0, 1);
for j=1:p
    
    fprintf(vertex_fid, '%d %1.16e %1.16e %d\n', j, x_points(j), y_points(j), 1);
end

fprintf(vertex_fid, '\n%d %d\n', p, 1);

for j=1:p

    fprintf(vertex_fid, '%d %d %d %d\n', j, facets(j,1), facets(j,2), 1);
    
end


fprintf(vertex_fid, '\n%d \n', 0);

%fprintf(vertex_fid, '%d %1.16e %1.16e %1.16e\n', 1, 0, 0, 0);



fclose(vertex_fid);


%
% spring_fid = fopen([mesh_name num2str(N) '.spring'], 'w');
% fprintf(spring_fid, '%d\n', springs_count);
% for s = 1:springs_count
%     fprintf(spring_fid, '%d %d %1.16e %1.16e\n', springs(s,1), springs(s,2), kappa_spring*ds/(ds^2), springs(s,3));
% end
% fclose(spring_fid);
%
%
%
% % beam_fid = fopen([mesh_name num2str(N) '.beam'], 'w');
% % fprintf(beam_fid, '%d\n', beams_count);
% % for b = 1: beams_count
% %     fprintf(beam_fid, '%d %d %d %1.16e\n', beams(b,1), beams(b,2), beams(b,3), kappa_beam*ds/(ds^4));
% % end
% % fclose(beam_fid);
%
% beam_fid = fopen([mesh_name num2str(N) '.beam'], 'w');
% fprintf(beam_fid, '%d\n', beams_count);
% for b = 1: beams_count
%     fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e %1.16e\n', beams(b,1), beams(b,2), beams(b,3), kappa_beam*ds/(ds^4), beams(b,4), beams(b,5), beams(b,6));
% end
% fclose(beam_fid);

%%

