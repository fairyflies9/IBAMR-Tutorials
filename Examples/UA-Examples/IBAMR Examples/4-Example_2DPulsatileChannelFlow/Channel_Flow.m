function Channel_Flow()

% Last Modified: 8/2/13
% Created By: Nick Battista (nickabattista@gmail.com)
%
% -This function constructs a channel geometry used to test time-dependent
%       and spatially dependent boundary conditions under the IBAMR 
%       framework.
%
% -It also prints out the necessary '.txt' files holding the vertex,
%       target, spring, and beam information.
% 
% -NOTE: It is important to match the geometry information (i.e.,
%        computational domain length L and d, with those in the 'input2d' file)


L = 5;          %Length of computational domain: [-L/2,L/2]
ds = L/(2*512); %Distance between Lagrangian pts.
d = 1;          %Channel Diameter

[X Y N_top] = get_Vertices(L,d,ds); %Constructs the channel geometry

plot_Geometry(L,X,Y) %Plots all the (X,Y) vertices

%Spring, Beam, Target Pt. Coefficients
kappa_target = 1e4;
kappa_spring = 1e-10;
kappa_beam = 1e-10;

%Actual Writing on Input Files for IBAMR
write_trab_vertices(X,Y);
write_trab_target(X,kappa_target,ds);
write_trab_spring(X,kappa_spring,ds,N_top);
write_trab_beam(X,Y,kappa_beam,ds,N_top)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: get_Vertices -> generates the geometry vertices.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X Y Ntop] = get_Vertices(L,d,ds)

%L = length of computational domain [-L/2,L/2]
%d = diamter of channel
%ds = spacing between Lagrangian points

x = -L/2 + ds/2; %To get first lagrangian point away from boundary of domain

i=0;
while x < L/2;
    i = i+1;
    X_top(i) = x;
    x = x + ds;
end

Ntop = i;

Y_top = d/2*ones(1,Ntop);

X_bot = X_top;

Y_bot = -Y_top;

X = [X_top X_bot]';
Y = [Y_top Y_bot]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: plot_Geometry -> plots all the (X,Y) vertices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Geometry(L,X,Y)

figure(1)
plot(X,Y,'*')
axis([-L/2 L/2 -L/2 L/2]);
title('Channel Flow Geometry');
xlabel('x');
ylabel('y');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write out VERTEX pts (x,y) for the trabeculae in straight tube
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_trab_vertices(x,y)

    N = length(x);

    vertex_fid = fopen(['channel' '.vertex'], 'w');

    fprintf(vertex_fid, '%d\n', N );

    %Prints Outer Tube Vertices
    for s = 0:N-1
        X_v = x(s+1);
        Y_v = y(s+1);
        fprintf(vertex_fid, '%1.16e %1.16e\n', X_v, Y_v);
    end

    fclose(vertex_fid); %%%ENDS HEART_TUBE VERTICES
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write out TARGET pts (x,y) for Lagrangian pts. in straight tube
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_trab_target(X,kappa_target,ds)

    N = length(X);
    
    %Write out the target point information
    target_fid = fopen(['channel' '.target'], 'w');

    fprintf(target_fid, '%d\n', N);  

    for s = 0:N-1
        fprintf(target_fid,'%d %1.16e\n',s,kappa_target*ds/(ds^2));
    end

    fclose(target_fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write out SPRING pts (x,y) for Lagrangian pts. in straight tube
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_trab_spring(X,kappa_spring,ds,N_top)

   N = length(X);   

   spring_fid = fopen(['channel' '.spring'], 'w');

   fprintf(spring_fid, '%d\n', N-2 ); 
   
   spring_force = kappa_spring*ds/(ds^2);

    %SPRINGS BETWEEN VERTICES
    for s = 0:N-2
        if s <= N_top-2
            fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, spring_force, ds);
        elseif s >= N_top
            fprintf(spring_fid, '%d %d %1.16e %1.16e\n', s, s+1, spring_force, ds);
        end
    end
    
    fclose(spring_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Write out BEAM pts (x,y) for Lagrangian pts. in straight tube
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_trab_beam(X,Y,kappa_beam,ds,N_t)

    N = length(X);
    
    beam_fid = fopen(['channel' '.beam'], 'w');

    fprintf(beam_fid, '%d\n', N-4 );

    %BEAMS BETWEEN OUTER TUBE VERTICES
    for s = 0:N-3     %Loops over all Vertices - 3
        
        C1= X(s+1) -2*X(s+2) + X(s+3); %Compute x-curvature
        C2= Y(s+1) -2*Y(s+2) + Y(s+3); %Compute y-curvature
        
        if s <= N_t-3    %Stops at 3rd to last vertex on top
            beam_force = kappa_beam*ds/(ds^4);
            fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', s, s+1, s+2, beam_force , C1, C2);
        elseif s >= N_t
            beam_force = kappa_beam*ds/(ds^4);
            fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', s, s+1, s+2, beam_force, C1, C2);
        end      
    end
    
    fclose(beam_fid);

