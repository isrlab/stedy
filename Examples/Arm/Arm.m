% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% This example demonstrates a tensegrity-based robotic arm with 16 bars and
% 18 strings. 
%
% We first generate associated parameters that are stored in tData, and we
% calculate an equilibrium force that is used in the simulation. An
% external force, which defined in 'ext_F.m', is applied to the structure.
% Plots and animation are made after the simulation.

clc; clear;

%% Define initial nodal coordinates
theta = pi/4; % To help describe angular position of any node
L = 1; % Length of all bars equal
vz = L*sin(theta); vx = L*cos(theta);
N = [0 0 0;
     0 0 -vz;
     vx 0 -vz;
     2*vx 0 -vz;
     2*vx 0 0;
     3*vx 0 -vz;
     4*vx 0 -vz;
     4*vx 0 0;
     5*vx 0 -vz;
     6*vx 0 0;
     5*vx 0 vz;
     4*vx 0 vz;
     3*vx 0 vz;
     2*vx 0 vz;
     vx 0 vz;
     0 0 vz;]';

nNodes = size(N,2); % Number of Nodes
nn = numel(N); % number of elements in N - 3*number of nodes: Probably not required



fixedNodes = [1 1 1;1 1 1; zeros(13,3);1 1 1]'; % To identify which coordinates of a node are fixed: 1-fixed, 0-unfixed
                                                % Same size as N

% Connectivity matrices
% C = [Cb; Cs;];
nBar = 16; % Number of bars
nStr = 18; % Number of strings
C = zeros(nBar+nStr,nNodes); % Connectivity Matrix C := [Cb; Cs];

% Presenting another way to define Cb and Cs
% Bars
C(1,1) = -1; C(1,3) = 1;
C(2,3) = -1; C(2,5) = 1;
C(3,5) = -1; C(3,4) = 1;
C(4,5) = -1; C(4,6) = 1;
C(5,6) = -1; C(5,8) = 1;
C(6,8) = -1; C(6,7) = 1;
C(7,8) = -1; C(7,9) = 1;
C(8,9) = -1; C(8,10) = 1;
C(9,10) = -1; C(9,11) = 1;
C(10,11) = -1; C(10,8) = 1;
C(11,8) = -1; C(11,12) = 1;
C(12,8) = -1; C(12,13) = 1;
C(13,13) = -1; C(13,5) = 1;
C(14,5) = -1; C(14,14) = 1;
C(15,5) = -1; C(15,15) = 1;
C(16,15) = -1; C(16,1) = 1;

% Strings
C(nBar+1,1) = -1; C(nBar+1,2) = 1;
C(nBar+2,2) = -1; C(nBar+2,3) = 1;
C(nBar+3,3) = -1; C(nBar+3,4) = 1;
C(nBar+4,4) = -1; C(nBar+4,6) = 1;
C(nBar+5,6) = -1; C(nBar+5,7) = 1;
C(nBar+6,7) = -1; C(nBar+6,9) = 1;
C(nBar+7,11) = -1; C(nBar+7,12) = 1;
C(nBar+8,12) = -1; C(nBar+8,13) = 1;
C(nBar+9,13) = -1; C(nBar+9,14) = 1;
C(nBar+10,14) = -1; C(nBar+10,15) = 1;
C(nBar+11,15) = -1; C(nBar+11,16) = 1;
C(nBar+12,16) = -1; C(nBar+12,1) = 1;
C(nBar+13,1) = -1; C(nBar+13,5) = 1;
C(nBar+14,5) = -1; C(nBar+14,8) = 1;
C(nBar+15,8) = -1; C(nBar+15,10) = 1;
C(nBar+16,9) = -1; C(nBar+16,11) = 1;
C(nBar+17,6) = -1; C(nBar+17,13) = 1;
C(nBar+18,3) = -1; C(nBar+18,15) = 1;

bars.index = 1:nBar; % To be used in constructing Cb
strings.index = nBar+1:nBar+nStr; % To be used in constructing Cs

Cb = C(bars.index,:); % Connectivity matrix for bars
Cs = C(strings.index,:); % Connectivity matrix for strings

tData = tensegConstruct(N,Cb,Cs,fixedNodes);

% Vector of Initial Nodal Positions and Velocities
x0 = [N(:);0*N(:)]; 

%% Material properties
% Strings
strings.r = 1/1000*ones(1,nStr);
strings.E = 2e9*ones(1,nStr); % Young's Modulus of Nylon
strings.rLP = 0.8*ones(1,nStr);  % Rest length 

% Bars
bars.rho = 1300*ones(1,nBar); 
bars.r = 0.01*ones(1,nBar);

% Point Masses
Mp = ones(1,tData.nPm); % All point masses initialised with a mass of 1

g = [0;0;-9.806]; % Gravity

tData = tensegGenMat(tData,bars,strings,Mp,g);

%% Simulation Inputs
tData.F = 1; % If 1, external force present in structure, else no external forces.
tData.Correction = 1; % If 1, constraint correction inclusive of total energy constraint If 0, only linear and bar length constraint violations corrected. 

tEnd = 20; % Simulation End Time

x0 = [x0;0]; % Initial Condition - [Position; Velocity; Work];

options = odeset('RelTol',1e-10,'AbsTol',1e-10); % ODE options

[t,y] = tensegSim(x0,tEnd,tData,options);

%% Equilibrium
% To find force densities in the bars and strings at equilibrium when
% under gravity and subjected to external forces. 
tData.minforce = 1; % Lower bound for force densities in the strings
tData = tensegEq(x0,N(:),tData);

%% Plotting

% Plot configuration
AZ = 0; % Azimuth angle in degrees
EL = 0; % Elevation angle in degrees
axLims = [-2.5 4.5 -4 4 -4 0]; 

plot_configuration(N(:),tData,AZ,EL,axLims);

% Plot Output Trajectories
plotMotion(t,y,tData);

% Plot Constraint Violations
plotConstr(t,y,tData);

% Animation
frameRate = 1/0.05;
time_loc = [1 0 2.5];
animateTenseg(t,y,tData,time_loc,'TEST_ARM.mp4','MPEG-4',frameRate,AZ,EL,axLims);

