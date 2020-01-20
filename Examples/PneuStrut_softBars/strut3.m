% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
% 
% This file demonstrates an example of a class-1 tensegrity structure with
% 3 bars and 9 strings.
%
% We first generate associated parameters that are stored in tData, and we
% calculate an equilibrium force that is used in the simulation. An
% external force, defined in 'ext_F.m', is applied to the structure.
% Plots and animation are made after the simulation.

clc; clear; 

%% Node generation
r32= sqrt(3)/2;
N = [0 -r32  r32  0 -r32 r32;
     1 -0.5 -0.5 -1  0.5 0.5;
     0    0    0  1    1   1;];

fixedNodes = zeros(size(N));
fixedNodes(:,1) = [1 1 1]';
fixedNodes(:,2) = [1 1 1]';
fixedNodes(:,3) = [1 1 1]';

% Connectivity Matrices
Cb = [-1  0  0 1 0 0;
       0 -1  0 0 0 1;
       0  0 -1 0 1 0;];

Cs = [-1  1  0  0  0  0;
       0 -1  1  0  0  0;
       1  0 -1  0  0  0;
       0  0  0 -1  1  0;
       0  0  0  0 -1  1;
       0  0  0  1  0 -1;
      -1  0  0  0  0  1;
       0 -1  0  0  1  0;
       0  0 -1  1  0  0;];

tData = tensegConstruct(N,Cb,Cs,fixedNodes);


% Vector of Initial Nodal Positions and Velocities
x0 = [N(:);0*N(:)];

%% Properties
% Strings
strings.index = tData.nBar+1:tData.nBar+tData.nStr;
strings.r = 1/1000*ones(1,tData.nStr);
strings.E = 2e9*ones(1,tData.nStr); % Young's Modulus of Nylon
strings.rLP = 1*ones(1,tData.nStr);  % Rest length percentage

% Bars
bars.index = 1:tData.nBar;
bars.rho = 960*ones(tData.nBar,1);
bars.r = 2/100*ones(tData.nBar,1);
bars.nu = 0.46*ones(tData.nBar,1); % Poisson's ratio (hdpe)
bars.E = 1E9*ones(tData.nBar,1); % Young's Modulus of hdpe
% bars.rLP = 0.8*ones(tData.nBar,1); % Rest Length of bars

% Point Masses
Mp = ones(1,tData.nPm);

g = [0;0;-9.806]; % Gravity

tData = tensegGenMat(tData,bars,strings,Mp,g);

%% Simulation Inputs
tData.F = 1; % If 1, external forces are present in structure, not if 0.
% tData.damper = 0.1*ones(1,tData.nStr);
% Dampers are not being used on strings in this example. 

%% Equilibrium
% To find force densities in the bars and strings at equilibrium when
% under gravity and subjected to external forces. 
tData.minforce = -Inf; % Lower bound for force densities in the strings
% Equilibrium for compressible bars, removed nonlinear constraints
tDataComp = tensegEqComp(x0,N(:),tData); 
tData = tensegEqComp(x0,N(:),tData);

%% Final Simulation
tEnd = 1; % Simulation End Time

options = odeset('RelTol',1e-10,'AbsTol',1e-10); % ODE options

x0 = [x0;0]; % Initial Condition - [Position; Velocity; Work];

[simTime,tInt] = tensegSimTime(options,tEnd);

tData.Correction = 2; % Rigid Bar 
tic
[t,y] = tensegSim(x0,simTime,tData,options);
compTime = toc

tDataComp.Correction = 3; % Compressible Bar 
tic
[tFlex,yFlex] = tensegSim(x0,simTime,tDataComp,options);
compTimeFlex = toc

%% Plotting

% Plot configuration
% AZ = 30;
% EL = 45;
% axLims = [min(N(1,:))-.3 max(N(1,:))+.3 min(N(2,:))-.3 max(N(2,:))+.3 min(N(3,:))-.3 max(N(3,:))+.3];
% plot_configuration(N(:),tData,AZ,EL,axLims);

% Plot Output Trajectories
% plotMotion(tFlex,yFlex,tDataComp);

% Diff. in Motion Plots
plotComp3strut_flex(t,y,tFlex,yFlex);

% print(figure(1),'comp3strutFlex_MotionNode4','-depsc');
% print(figure(2),'comp3strutFlex_MotionNode5','-depsc');
% print(figure(3),'comp3strutFlex_MotionNode6','-depsc');

% Plot Constraint Violations
% plotConstr(t,y,tData);

% Animation
% filename = 'Ball_Animation';
% formatSpec = 'MPEG-4';
% frameRate = 1/tInt;
% time_loc = [0 1.8 0];
% animateTenseg(t,y,tData,time_loc,filename,formatSpec,frameRate,AZ,EL,axLims);