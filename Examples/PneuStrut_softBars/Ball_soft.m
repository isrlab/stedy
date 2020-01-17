% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
% 
% This file demonstrates an example of a ball-shaped tensegrity strucute,
% which is formulated by 6 bars, 32 strings and 1 point mass at the geometric 
% center of the structure.
%
% We first generate associated parameters that are stored in tData, and we
% calculate an equilibrium force that is used in the simulation. An
% external force, defined in 'ext_F.m', is applied to the structure.
% Plots and animation are made after the simulation.

clc; clear; 

%% Node generation
% First, we need to have the coordinates of a regular icosahedron, which is
% formulated from the golden ratio, and the coordinates are described by
% (0,+-1,+-phi), where phi is the golden ratio. The payload is then at the
% origin.  In addition, we would like to have 3 nodes (node 1,2,9) laying
% on x-y plane. To do this, we move the origin to node 1 and rotate the
% body around node 1 in angle alpha such that 3 nodes are on the plane.

phi = (1+sqrt(5))/2; % Golden ratio
beta = 6; % There are 6 bars.
temp = [0 -1 phi;0 -1 -phi;0 1 phi;0 1 -phi;]/(2*phi); % Change the bar length to 1 meter; 
N = [temp;temp(:,[3 1 2]);temp(:,[2 3 1]);0 0 0]'; % Original order (too complicated and hard to figure out)
N = N(:,[2 4 6 8 10 12 1 3 5 7 9 11 13]); % Rearrange the order such that first beta nodes are the bases and last beta nodes are the terminal nodes.
N = N-kron(ones(1,size(N,2)),N(:,1)); % Move the origin to node 1.
fun = @(alpha)(sin(alpha)*N(1,9)+cos(alpha)*N(3,9))*100000;
alpha=fsolve(fun,0,optimoptions('fsolve','Display','off')); % Calculate alpha such that the position in z of node 9 is zero.
dcm = [cos(alpha) 0 -sin(alpha);0 1 0;sin(alpha) 0 cos(alpha)]; % Direction cosine matrix
N = dcm*N; % Rotate the coordinate about the original so that node 1 2 9 is on the ground (z=0).
Nd = zeros(size(N));

fixedNodes = zeros(size(N));
fixedNodes(:,1) = [1 1 1]';
fixedNodes(:,2) = [1 1 1]';
fixedNodes(:,9) = [1 1 1]';

% Connectivity Matrices
Cb = [-eye(beta) eye(beta) zeros(beta,1)];

Cs = [-1     0     1     0     0     0     0     0     0     0     0     0     0;
      -1     0     0     0     1     0     0     0     0     0     0     0     0;
      -1     0     0     0     0     1     0     0     0     0     0     0     0;
      -1     0     0     0     0     0     0     0     1     0     0     0     0;
       0    -1     1     0     0     0     0     0     0     0     0     0     0;
       0    -1     0     0     0     0     0     0     1     0     0     0     0;
       0    -1     0     0     0     0     0     0     0     0     1     0     0;
       0    -1     0     0     0     0     0     0     0     0     0     1     0;
       0     0    -1     0     1     0     0     0     0     0     0     0     0;
       0     0    -1     0     0     0     0     0     0     0     1     0     0;
       0     0     0    -1     1     0     0     0     0     0     0     0     0;
       0     0     0    -1     0     0     1     0     0     0     0     0     0;
       0     0     0    -1     0     0     0     1     0     0     0     0     0;
       0     0     0    -1     0     0     0     0     0     0     1     0     0;
       0     0     0     0    -1     0     1     0     0     0     0     0     0;
       0     0     0     0     0    -1     1     0     0     0     0     0     0;
       0     0     0     0     0    -1     0     0     1     0     0     0     0;
       0     0     0     0     0    -1     0     0     0     1     0     0     0;
       0     0     0     0     0     0    -1     0     0     1     0     0     0;
       0     0     0     0     0     0     0    -1     0     1     0     0     0;
       0     0     0     0     0     0     0    -1     0     0     1     0     0;
       0     0     0     0     0     0     0    -1     0     0     0     1     0;
       0     0     0     0     0     0     0     0    -1     0     0     1     0;
       0     0     0     0     0     0     0     0     0    -1     0     1     0;
       0     0     1     0     0     0     0     0     0     0     0     0    -1;
       0     0     0     1     0     0     0     0     0     0     0     0    -1;
       0     0     0     0     1     0     0     0     0     0     0     0    -1;
       0     0     0     0     0     1     0     0     0     0     0     0    -1;
       0     0     0     0     0     0     0     0     1     0     0     0    -1;
       0     0     0     0     0     0     0     0     0     1     0     0    -1;
       0     0     0     0     0     0     0     0     0     0     1     0    -1;
       0     0     0     0     0     0     0     0     0     0     0     1    -1];

tData = tensegConstruct(N,Cb,Cs,fixedNodes);


% Vector of Initial Nodal Positions and Velocities
x0 = [N(:);Nd(:)];

%% Properties
% Strings
strings.index = tData.nBar+1:tData.nBar+tData.nStr;
strings.r = 1/1000*ones(1,tData.nStr);
strings.E = 2e9*ones(1,tData.nStr); % Young's Modulus of Nylon
strings.rLP = 0.9*ones(1,tData.nStr);  % Rest length percentage

% Bars
bars.index = 1:tData.nBar;
bars.rho = 960*ones(tData.nBar,1);
bars.r = 2/100*ones(tData.nBar,1);
bars.nu = 0.46*ones(tData.nBar,1); % Poisson's ratio (HDPE)
bars.E = 1E9*ones(tData.nBar,1); % Young's Modulus of HDPE

% Point Masses
Mp = ones(1,tData.nPm);

g = [0;0;-9.806]; % Gravity

tData = tensegGenMat(tData,bars,strings,Mp,g);

%% Simulation Inputs
tData.F = 1; % If 1, external forces are present in structure, not if 0.
tData.Correction = 3; % If 1, constraint correction inclusive of total energy constraint If 0, only linear and bar length constraint violations corrected. 
% tData.damper = 0.1*ones(1,tData.nStr);
% Dampers are not being used on strings in this example. 

%% Equilibrium
% To find force densities in the bars and strings at equilibrium when
% under gravity and subjected to external forces. 
tData.minforce = -Inf; % Lower bound for force densities in the strings
tData = tensegEqComp(x0,N(:),tData);

%% Final Simulation
tEnd = 1; % Simulation End Time

options = odeset('RelTol',1e-10,'AbsTol',1e-10); % ODE options

x0 = [x0;0]; % Initial Condition - [Position; Velocity; Work];

[simTime,tInt] = tensegSimTime(options,tEnd);

[t,y] = tensegSim(x0,simTime,tData,options);

%% Plotting

% Plot configuration
AZ = 30;
EL = 45;
axLims = [min(N(1,:))-.3 max(N(1,:))+.3 min(N(2,:))-.3 max(N(2,:))+.3 min(N(3,:))-.3 max(N(3,:))+.3];
plot_configuration(N(:),tData,AZ,EL,axLims);

% Plot Output Trajectories
plotMotion(t,y,tData);

% Plot Constraint Violations
% plotConstr(t,y,tData);

% Animation
% filename = 'Ball_Animation';
% formatSpec = 'MPEG-4';
% frameRate = 1/tInt;
% time_loc = [0 1.8 0];
% animateTenseg(t,y,tData,time_loc,filename,formatSpec,frameRate,AZ,EL,axLims);