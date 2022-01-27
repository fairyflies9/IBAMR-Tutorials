%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Test Geometry/Spring Connections for IB2d/IBAMR geometry
%       
%      Purpose: (1) Plots Lagrangian points themselves
%               (2) Plots Lagrangian points & spring connections
%
%      NOTE: user will need to specify:
%
%               [1] Line 23: whether it is IBAMR input files since indexing 
%                   starts  at 0 in C++ (note: pyIB2d can also start at 0.)
%
%               [2] Line 40: How many columns are in the .spring file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_IB2d_IBAMR_Geometry()


%------------------------------------
% FLAG FOR IB2d or IBAMR
%------------------------------------
flag_IBAMR = 0; % 0 for IB2d, 1 for IBAMR


%------------------------------------
% Get xy-VERTICES:
%     --> XY: 1st column x-Pts
%             2nd column y-Pts
%------------------------------------
fileName = 'rubberband.vertex';
XY = read_Vertex_Data(fileName);

%------------------------------------------------------------
% Get spring connections:
%     --> springs: 1st column: leader lag. pt. index
%                  2nd column: follower lag. pt. index
%------------------------------------------------------------
fileName = 'rubberband.spring';   % name of .spring file
numCols = 4;                      % number of columns in .spring file
springs = read_Spring_Data(fileName,numCols);

%**********************************************************
% NOTE IF USING IBAMR ADD +1 to COL 1 & COL 2 OF SPRINGS
%       FOR C++ INDEXING 
%**********************************************************
if flag_IBAMR == 1
    ind_Lags = 0:length(XY(:,1))-1;      % Indices of lag pts (starts at 0)
    springs(:,1:2) = springs(:,1:2) + 1; % Correct indices of spring connections
else
    ind_Lags = 1:length(XY(:,1));        % Indices of lag pts (starts at 1)
end


%--------------------------------------------------------------
% Explicitly define leader and follwer indices for plotting
%--------------------------------------------------------------
inds_L = springs(:,1); % Leader indices of springs
inds_F = springs(:,2); % Follower indices of springs


%---------------------------------
% FIG 1: Test Vertices
%---------------------------------
figure(1)
for i=1:length(XY(:,1))
    plot(XY(i,1),XY(i,2),'k.','MarkerSize',20); hold on;
    text(XY(i,1),XY(i,2),['  ' num2str(ind_Lags(i))] );
end
xlabel('x');
ylabel('y');
set(gca,'FontSize',22);

%-------------------------------------------------
% FIG 2: Test Vertices + Spring Connections
%-------------------------------------------------
figure(2)
plot(XY(:,1),XY(:,2),'k.','MarkerSize',20); hold on;
%
for i=1:length(inds_L)
    
    % Starting and ending x-Values of i-th spring
    xSpr = [XY(inds_L(i),1) XY(inds_F(i),1)];
    
    % Starting and ending y-Values of i-th spring
    ySpr = [XY(inds_L(i),2) XY(inds_F(i),2)];
    
    % Plot the i-th spring (either red of blue, just to change up
    %                       successive coloring for sole purpose of testing)
    if mod(i,2)==1
        plot(xSpr,ySpr,'r-','LineWidth',2); hold on;
    else
        plot(xSpr,ySpr,'b-','LineWidth',2); hold on;
    end
    
end
%
for i=1:length(XY(:,1))
    plot(XY(i,1),XY(i,2),'k.','MarkerSize',20); hold on;
    text(XY(i,1),XY(i,2),['  ' num2str(ind_Lags(i))] );
end
%
xlabel('x');
ylabel('y');
set(gca,'FontSize',22);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of vertex pts and all the vertex pts from the
%           .vertex file 
%
%   NOTE: pass string of input file name, e.g., 'rubberband.vertex' 
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DATA = read_Vertex_Data(filename)

fileID = fopen(filename);

% Read in the file, use 'CollectOutput' to gather all similar data together
% and 'CommentStyle' to to end and be able to skip lines in file.
C = textscan(fileID,'%f %f','CollectOutput',1);

fclose(fileID);     %Close the data file.

vertices = C{1};    %Stores all read in data in vertices (N+1,2) array

N = vertices(1,1);  % # of Lagrangian Pts 
DATA = vertices(2:end,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the # of springs and all of the spring information
%               from the .spring file
%
%   NOTE: pass string of input file name, e.g., 'rubberband.spring' 
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DATA = read_Spring_Data(filename,numCols)

fileID = fopen(filename);

vecRead = '%f';
for i=1:numCols-1
   vecRead = [vecRead ' %f']; 
end

% Read in the file, use 'CollectOutput' to gather all similar data together
% and 'CommentStyle' to to end and be able to skip lines in file.
C = textscan(fileID,vecRead,'CollectOutput',1);

fclose(fileID);     %Close the data file.

vertices = C{1};    %Stores all read in data in vertices (N+1,2) array

N = vertices(1,1);  % # of Springs
DATA = vertices(2:end,:);

