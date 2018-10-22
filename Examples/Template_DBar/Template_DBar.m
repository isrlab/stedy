% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

clc; clear;

% This is a template file to help the user understand how
% inputs are to be provided to the functions inside the software.

% The structure being simulated is a D-bar. 

%% Define Initial Nodal Coordinates and Connectivity Matrices
N = [0      1     1     0;
     0      0     0     0;
     0      0     1     1;];
% Structure of N:
% N = [N1  N2  N3 ..]
%      N1x N2x N3x ..;  
%      N1y N2y N3y ..;
%      N1z N2z N3z ..;];

fixedNodes = [1 1 1;zeros(3,3)]'; 
% To identify which coordinates of a node are fixed
% 1 - fixed, 0 - free
% Matrix of the same size as N

% Connectivity Matrices
C = [ 1  0  0 -1;
     -1  1  0  0;
      0 -1  1  0;
      0  0 -1  1;
      0 -1  0  1;
     -1  0  1  0;]; % Connectivity Matrix C := [Cb; Cs];
nBar = 4; % Number of bars
nStr = size(C,1) - nBar; % Number of strings

bars.index = 1:nBar; % To be used in constructing Cb
strings.index = nBar+1:nBar+nStr; % To be used in constructing Cs

Cb = C(bars.index,:); % Connectivity matrix for bars
Cs = C(strings.index,:); % Connectivity matrix for strings

% It is not necessary to define the connectivity matrices Cb, Cs this way.
% It is perfectly acceptable to define Cb and Cs independently without
% constructing C as C = [Cb;Cs;]; 
% Point masses are placed exclusively at those nodes that do not connect to
% any bars. 

tData = tensegConstruct(N,Cb,Cs,fixedNodes);
% Function to create a tensegrity structure with the prescribed data.
% tData shall be the MATLAB structure containing all the important fields
% that facilitate the dynamics computations.  

x0 = [N(:);0*N(:)];
% x0 is the vector of initial nodal coordinates on top of the vector of
% initial velocities. Velocities have been chosen to be 0 in all the given
% examples. 

%% Material and Simulation Environment Properties
% Strings
strings.r = 1/1000*ones(1,tData.nStr); % Radius of strings
strings.E = 7e7*ones(1,tData.nStr); % Young's Modulus of Nylon
strings.rLP = [1 0.7];  % Rest length of the strings
% 0.9 rest length means 90% of the current length of the string is its rest
% length. 
% Strings are modeled as springs but unlike springs, they can only exert
% tensile force, not compressive. 

% Bars
bars.rho = 500*ones(1,tData.nBar); % Density of bars
bars.r = 1/100*ones(1,tData.nBar); % Radius of bars. 
% Bars are modeled as cylinders with a specified density and radius. 

% Point Masses
Mp = 1/100*ones(1,tData.nPm); 
% Any and all point masses have been initialised with a mass of 1.
% No point masses in this structure. 

% Gravity
g = [0;0;0;];
% The gravity vector can be described any way the user chooses. In the
% T-bar example, gravity is taken to be zero. 

tData = tensegGenMat(tData,bars,strings,Mp,g);

%% Simulation Inputs
tData.F = 0; 
% If 1, external force present in structure, else no external forces.

tData.Correction = 1; 
% If 1, constraint correction inclusive of total energy constraint.
% If 0, only linear and bar length constraint violations corrected. 

tData.damper = ones(1,tData.nStr)*0.1; 
% All strings initialised with a damper whose damping coefficient is 1. 
% The user is free to choose which strings to damp by assigning an 
% appropriate coefficient to the corresponding element in the vector. 

%% Equilibrium Calculations
% tData.minforce = eps; % The minimum force density (scalar) in the equilibrium position.
% tData = tensegEq(x0,N(:),tData);

% It is not requried to find the equilibrium force densities in this
% structure, and hence, it is not being run in this simulation. For
% examples, consider Test_Ball.m and Test_Arm.m. It is required of the user
% to define external forces in ext_F.m and extF_eq.m accordingly. 

%% Final Simulation
x0 = [x0;0;];
% The last element in the vector describes the work done on the structure
% at t=0. 

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
% ODE options. The user is encouraged to read MATLAB's ODE45 documentation
% to learn how to provide these inputs. 

tEnd = 10; % Simulation End Time

[simTime,tInt] = tensegSimTime(options,tEnd);

[t,y] = tensegSim(x0,simTime,tData,options);

%% Plotting 

%% Plot Initial Configuration
AZ = 0; % Azimuth angle in degrees
EL = 0; % Elevation angle in degrees
axLims = [-5 5 -3 3 -5 5]; % Axis limits in the figure window for animateTenseg
plot_configuration(N(:),tData,AZ,EL,axLims);
% The user is recommended to look at MATLAB documentation to understand how
% to provide AZ and EL inputs. The final input in the above command
% represents the dimensions of the figure window in which the plot is
% drawn. 

%% Trajectories
plotMotion(t,y,tData);
% Function to plot the trajectories of all the nodes in the tensegrity
% structure. 

%% Constraint Violations
plotConstr(t,y,tData);
% Function to plot both the bar length constraint violations as well as
% errors in total energy conservation generated by numerical integration. 

%% Animation
% AZ, EL and axLims shall be the same as used in plot_configuration
filename = 'DBar_Animation';
formatSpec = 'MPEG-4';
frameRate = 1/tInt; % Video Framerate
time_loc = [0 0 4]; % Location of displayed time, can be chosen differently by the user
animateTenseg(t,y,tData,time_loc,filename,formatSpec,frameRate,AZ, EL,axLims);
% The user is encouraged to look further at MATLAB documentation on
% VideoWriter to see what video formats are allowed. 

%% Saving Plots
% print(figure(1),'DBAR_Configuration','-djpeg');
% print(figure(2),'DBAR_MotionNode1','-djpeg');
% print(figure(3),'DBAR_MotionNode2','-djpeg');
% print(figure(4),'DBAR_MotionNode3','-djpeg');
% print(figure(5),'DBAR_MotionNode4','-djpeg');
% print(figure(6),'DBAR_Bar1LengthError','-djpeg');
% print(figure(7),'DBAR_Bar2LengthError','-djpeg');
% print(figure(8),'DBAR_Bar3LengthError','-djpeg');
% print(figure(9),'DBAR_Bar4LengthError','-djpeg');
% print(figure(10),'DBAR_EnergyError','-djpeg');
