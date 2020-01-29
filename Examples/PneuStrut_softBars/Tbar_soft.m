% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
% 
% To compare Tbar w/o correction against Tbar_flex.
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
bars.rho = 960*ones(1,tData.nBar); % Density of bars
bars.nu = 0.46*ones(1,tData.nBar); % Poisson's ratio of bars (HDPE)
bars.E = 1e9*ones(1,tData.nBar); % Young's modulus of bars (HDPE)

% Point Masses
Mp = ones(1,tData.nPm); % All point masses initialised with a mass of 1

g = [0;0;0]; % Gravity

tData = tensegGenMat(tData,bars,strings,Mp,g);

% The structure is in equilibrium in the initial position. Hence, there
% is no need to call the function tensegEq.

%% Simulation Inputs
tData.F = 0; % If 1, external forces are present in structure, not if 0.
tData.damper = ones(1,tData.nStr); % All strings initialised with dampers whose damping coefficient is 1. 

%% Final Simulation
tEnd = 10; % Simulation End Time

x0 = [x0;0]; % Initial Condition - [Position; Velocity; Energy];

options = odeset('RelTol',1e-10,'AbsTol',1e-10,'Refine',1);

[simTime,tInt] = tensegSimTime(options,tEnd);

tData.Correction = 3; % Compressible Bar 
tic
linSys = linSysComp(x0, tData);
% [tFlex,yFlex] = tensegSim(x0,simTime,tData,options);
compTimeFlex = toc
%% Plotting 
% Plot Output Trajectories
% plotMotion(tFlex, yFlex, tData);

% Overlay Motion Plots
% plotCompTbar_flex(t,y,tFlex,yFlex,tData);

% print(figure(3),'compTBarFlex_MotionNode3','-depsc');
% print(figure(4),'compTBarFlex_MotionNode4','-depsc');

