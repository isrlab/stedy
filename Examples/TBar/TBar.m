% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

clc; clear;

%% Define initial nodal coordinates and connectivity matrices
theta = pi/4; % To help describe angular position of any node

N = [0  5*cos(theta)  5*cos(theta)              0;
     0            0              0              0;
     0            0   5*sin(theta)   5*sin(theta);];

fixedNodes = [1 1 1;1 1 1; 0 1 0; 0 1 0;]'; % To identify which coordinates of a node are fixed: 1-fixed, 0-unfixed
                                            % Same size as N

% Connectivity Matrices
% Cb an Cs have been defined independently. 
Cb = [-1  0  1  0; 
       0 -1  0  1;];
  
Cs = [-1  1  0  0;
       0 -1  1  0;
       0  0 -1  1;
       1  0  0 -1;]; 

tData = tensegConstruct(N,Cb,Cs,fixedNodes);

% Vector of Initial Nodal Positions and Velocities
x0 = [N(:); 0*N(:)];

%% Material and Simulation Environment properties
% Strings
strings.r = 1/1000*ones(1,tData.nStr); % Radius of strings
strings.E = 7e7*ones(1,tData.nStr); % Young's Modulus of Nylon
strings.rLP = [1 0.7 1 0.7];  % Rest lengths of the strings: 0.7 means 70% 

% Bars
bars.r = 0.05*ones(1,tData.nBar); % Radius of bars
bars.rho = 500*ones(1,tData.nBar); % Density of bars

% Point Masses
Mp = ones(1,tData.nPm); % All point masses initialised with a mass of 1

g = [0;0;0]; % Gravity

tData = tensegGenMat(tData,bars,strings,Mp,g);

% The structure is in equilibrium in the initial position. Hence, there
% is no need to call the function tensegEq.

%% Simulation Inputs
tData.F = 0; % If 1, external forces are present in structure, not if 0.
tData.Correction = 1; % If 1, constraint correction inclusive of total energy constraint If 0, only linear and bar length constraint violations corrected. 
tData.damper = ones(1,tData.nStr); % All strings initialised with dampers whose damping coefficient is 1. 

%% Final Simulation
tEnd = 10; % Simulation End Time

x0 = [x0;0]; % Initial Condition - [Position; Velocity; Energy];

options = odeset('RelTol',1e-10,'AbsTol',1e-10,'Refine',1);

[simTime,tInt] = tensegSimTime(options,tEnd);

[t,y] = tensegSim(x0,simTime,tData,options);

%% Plotting 

% Plot Configuration
AZ = 0; % Azimuth angle in degrees
EL = 0; % Elevation angle in degrees
axLims = [-6 6 -6 6 -8 8]; % Axis Limits in the figure window

plot_configuration(N(:),tData,AZ,EL,axLims);

% Plot Output Trajectories
plotMotion(t,y,tData);

% Plot Constraint Violations
plotConstr(t,y,tData);

% Animation
filename = 'TBar_Animation';
formatSpec = 'MPEG-4';
frameRate = 1/tInt;
time_loc = [0.5 0 5];
animateTenseg(t,y,tData,time_loc,filename,formatSpec,frameRate,AZ,EL,axLims);
