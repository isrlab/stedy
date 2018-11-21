% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
% 
% This file provides a test case for the user to check their version of the
% dynamics simulation with the results of a minimum-coordinates
% realisation, conventionally accepted as the ground truth. 
% 
clc; clear;

%% Define initial nodal coordinates and connectivity matrices
theta = pi/4; % To help describe angular position of any node
lb = 5;

N = [0  lb*sin(theta);
     0            0;  
     0  -lb*cos(theta)];

fixedNodes = [1 1 1;0 1 0;]'; % To identify which coordinates of a node are fixed: 1-fixed, 0-unfixed
                                            % Same size as N

% Connectivity Matrices
% Cb an Cs have been defined independently. 
Cb = [-1  1];
  
Cs = [];

tData = tensegConstruct(N,Cb,Cs,fixedNodes);

% Vector of Initial Nodal Positions and Velocities
x0 = [N(:); 0*N(:)];

%% Material and Simulation Environment properties
% Strings
strings.r = 1/1000*ones(1,tData.nStr); % Radius of strings
strings.E = 7e7*ones(1,tData.nStr); % Young's Modulus of Nylon
strings.rLP = [];  % Rest lengths of the strings: 0.7 means 70% 

% Bars
bars.r = 0.05*ones(1,tData.nBar); % Radius of bars
bars.rho = 500*ones(1,tData.nBar); % Density of bars

% Point Masses
Mp = ones(1,tData.nPm); % All point masses initialised with a mass of 1

g = [0;0;-9.806]; % Gravity

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

options = odeset('RelTol',1e-10,'AbsTol',1e-10);

tInt = 0.01;
simTime = 0:tInt:tEnd;

[t,y] = tensegSim(x0,simTime,tData,options);

%% MinReal
q_in = [theta;0]; % Initial Conditions

m = bars.rho*pi*bars.r^2*lb;
I = m*lb^2/3;

simTime = 0:tInt:tEnd;

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,q] = ode45(@Pendulum_MinRealODE,simTime,q_in,options,tData);

X = lb*sin(q(:,1)); Z = -lb*cos(q(:,1));

%% Comparing Results

% Plotting Comparison with Minimum Realization
figure(1); clf;

subplot(2,1,1);
plot(t,X-y(:,4)); hold on;
xlabel('Time (s)'); ylabel('x(t)'); 
set(gca,'box','off');

subplot(2,1,2);
plot(t,Z - y(:,6)); hold on;
xlabel('Time (s)'); ylabel('z(t)'); 
set(gca,'box','off');

% Test Constraint Violations
testConstraint(t,y,tData);
%% Display Results
if(norm(X-y(:,4)) && norm(Z-y(:,6)) < 1e-6)
    fprintf('Motion Error under 1e-6.\n');
else
    fprintf('Motion Error over 1e-6. Test failed.\n');
end

