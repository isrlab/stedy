% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
% 
% Plate (3 Tbars arranged in series) for linear systems' analysis.
clc; clear;

%% Define initial nodal coordinates and connectivity matrices
theta = pi/4; % To help describe angular position of any node

N1 = [0  5*cos(theta)  5*cos(theta)              0; 
      0            0              0              0;
      0            0   5*sin(theta)   5*sin(theta);];
N2 = [10*cos(theta) 15*cos(theta) 15*cos(theta) 10*cos(theta);
                  0             0             0             0;
                  0             0  5*sin(theta)  5*sin(theta);];
N = [N1 N2];              
fixedNodes1 = [1 1 0;1 1 0; 0 1 0; 0 1 0;]; % To identify which coordinates of a node are fixed: 1-fixed, 0-unfixed
                                            % Same size as N
fixedNodes2 = fixedNodes1;
fixedNodes = [fixedNodes1; fixedNodes2]';

% Connectivity Matrices
% Cb an Cs have been defined independently. 
Cb = zeros(6,8);
Cb(1,1) = -1; Cb(1,3) = 1;
Cb(2,2) = -1; Cb(2,4) = 1;
Cb(3,2) = -1; Cb(3,8) = 1;
Cb(4,5) = -1; Cb(4,3) = 1;
Cb(5,5) = -1; Cb(5,7) = 1;
Cb(6,6) = -1; Cb(6,8) = 1;
  
Cs = zeros(10,8);
Cs(1,1) = -1; Cs(1,2) = 1;
Cs(2,2) = -1; Cs(2,3) = 1;
Cs(3,3) = -1; Cs(3,4) = 1;
Cs(4,4) = -1; Cs(4,1) = 1;
Cs(5,2) = -1; Cs(5,5) = 1;
Cs(6,5) = -1; Cs(6,6) = 1;
Cs(7,6) = -1; Cs(7,7) = 1;
Cs(8,7) = -1; Cs(8,8) = 1;
Cs(9,3) = -1; Cs(9,8) = 1;
Cs(10,5) = -1; Cs(10,8) = 1;


tData = tensegConstruct(N,Cb,Cs,fixedNodes);

% Vector of Initial Nodal Positions and Velocities
v0 = [0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0;
      1 1 1 1 1 1 1 1;];
x0 = [N(:); v0(:);];
x0 = [N(:);0*N(:)];

%% Material and Simulation Environment properties
% Strings
strings.r = 1/1000*ones(1,tData.nStr); % Radius of strings
strings.E = 7e7*ones(1,tData.nStr); % Young's Modulus of Nylon
strings.rLP = ones(1,tData.nStr);  % Rest lengths of the strings: 0.7 means 70% 

% % Bars (Soft)
bars.r = 0.05*ones(1,tData.nBar); % Radius of bars
bars.rho = 960*ones(1,tData.nBar); % Density of bars
bars.nu = 0.46*ones(1,tData.nBar); % Poisson's ratio of bars (HDPE)
bars.E = 1e9*ones(1,tData.nBar); % Young's modulus of bars (HDPE)

% Bars (Metallic)
bars.r(end) = 0.05; % Radius of bars
bars.rho(end) = 2700; % Density of bars
bars.nu(end) = 0.30; % Poisson's ratio of bars (aluminium)
bars.E(end) = 200e9; % Young's modulus of bars

% Point Masses
Mp = ones(1,tData.nPm); % All point masses initialised with a mass of 1

g = [0;0;0]; % Gravity

tData = tensegGenMat(tData,bars,strings,Mp,g);

% The structure is in equilibrium in the initial position. Hence, there
% is no need to call the function tensegEq.

%% Simulation Inputs
tData.F = 0; % If 1, external forces are present in structure, not if 0.
% tData.damper = ones(1,tData.nStr); % All strings initialised with dampers whose damping coefficient is 1. 
tData.minforce = 1; % Lower bound for force densities in the strings
tData = tensegEqComp(x0,tData);

%% Final Simulation
% tEnd = 10; % Simulation End Time

x0 = [x0;0]; % Initial Condition - [Position; Velocity; Energy];

options = odeset('RelTol',1e-10,'AbsTol',1e-10,'Refine',1);

% [simTime,tInt] = tensegSimTime(options,tEnd);

tData.Correction = 3; % Compressible Bar 
% tic
[linSys,sysSS] = linSysCompDescriptor(x0, tData);
Asys = linSys.A; Bsys = linSys.Bu; 
Csys = linSys.C; Dsys = linSys.D;
% [tFlex,yFlex] = tensegSim(x0,simTime,tData,options);
minSys = minreal(sysSS);
Amin = minSys.A;
szAmin = size(Amin)
rMin = rank(Amin)
e = eig(Amin)
maxE = max(abs(e))
nnzE = nnz(abs(e)>1e-5)

%% Plotting 
% Plot Configuration
% AZ = 0; % Azimuth angle in degrees
% EL = 0; % Elevation angle in degrees
% axLims = [-20 20 -6 6 -20 20]; % Axis Limits in the figure window
% 
% plot_configuration(N(:),tData,AZ,EL,axLims);
% Plot Output Trajectories
% plotMotion(tFlex, yFlex, tData);

% Overlay Motion Plots
% plotCompTbar_flex(t,y,tFlex,yFlex,tData);

% print(figure(3),'compTBarFlex_MotionNode3','-depsc');
% print(figure(4),'compTBarFlex_MotionNode4','-depsc');

