%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nickabattista@gmail.com
% Date Created: May 27th, 2015
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "Hill+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nickabattista@gmail.com) know.
%
%--------------------------------------------------------------------------------------------------------------------%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: updates the target point positions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function targets = update_Target_Point_Positions(dt,current_time,targets)


IDs = targets(:,1);                 % Stores Lag-Pt IDs in col vector
%%%%%%%%%%%These are the most recent x- and y-coordinates for the target
%%%%%%%%%%%points. Same for the values of kStiffs
xPts= targets(:,2);                 % Previous x-Values of x-Target Pts.
yPts= targets(:,3);                  % Previous y-Values of y-Target Pts.
kStiffs = targets(:,4);             % Stores Target Stiffnesses 
N_target = length(targets(:,1));    % Gives total number of target pts!

%
% PARAMETERS
amp= 1/8;      %amplitude for the pump
freq= 1;     %frequency of the pump
stop_time = 6;  %stop actively moving beams at this time
first_bottom = 1;   %first target point on bottom beam
last_bottom = 64;   %last target point on bottom beam
first_top = 257;   %first target point on top beam
last_top = 320;   %last target point on top beam
%
% READ IN ORIGINAL yPT and xPT POSITIONS
%
yPTS = read_In_yPT_Positions('two_beams_512.vertex');
xPTS = read_In_xPT_Positions('two_beams_512.vertex');
%
% Create the pumping behavior
% To change the target point position, we will modify the target array. The
% y-position is the third column.

while current_time<stop_time
    for i=first_bottom:last_bottom
      targets(i,3) = yPTS(i) - amp*sin( 2*pi*freq*current_time );
    end
    for i=first_top:last_top
      targets(i,3) = yPTS(i) + amp*sin( 2*pi*freq*current_time );
    end
end
%yPts(N_top+1:end) = yPts(N_top+1:end) + 0.95*(G/2)*sin( 2*pi*freq*current_time );





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of vertex pts and all the vertex pts from the
%           .vertex file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PTS = read_In_yPT_Positions(struct_name)


filename = struct_name;  %Name of file to read in
fileID = fopen(filename);

% Read in the file, use 'CollectOutput' to gather all similar data together
% and 'CommentStyle' to to end and be able to skip lines in file.
C = textscan(fileID,'%f %f','CollectOutput',1);

fclose(fileID);     %Close the data file.

vertices = C{1};    %Stores all read in data in vertices (N+1,2) array

PTS = vertices(2:end,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PTS = read_In_xPT_Positions(struct_name)


filename = struct_name;  %Name of file to read in
fileID = fopen(filename);

% Read in the file, use 'CollectOutput' to gather all similar data together
% and 'CommentStyle' to to end and be able to skip lines in file.
C = textscan(fileID,'%f %f','CollectOutput',1);

fclose(fileID);     %Close the data file.

vertices = C{1};    %Stores all read in data in vertices (N+1,2) array

PTS = vertices(2:end,1);
