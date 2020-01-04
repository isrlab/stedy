function [t,y] = tensegSim(x0,simTime,tData,options)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
% 
% Function to simulate the tensegrity structure's dynamics in the
% prescribed environment
%
% INPUTS: [x0, tEnd, tData, options]
% 
% x0: Vector at t=0 := [Nodal coordinates; Nodal velocites; Work done];
%
% tEnd: Simulation End Time
%
% options: ODE45 options. The user is encouraged to look at MATLAB's ode45
% documentation to understand how to provide this input.
%
% NOTE: If 'Refine' option is not turned on in the ode options provided,
% the software will prompt the user to enter a suitable output interval
% time.
% 
% OUTPUTS: [t,y]
%
% t: Evaluation points, returned as a column vector. 
%
% y: Solutions returned as a matrix. Each row in the solution matrix y 
% corresponds to values returned in column vector t. 


nn = numel(tData.N);

if(tData.Correction == 0 || tData.Correction == 1) % Constraint Correction turned on
    tData.Energy = 0.5*x0(nn+1:end-1)'*tData.M*x0(nn+1:end-1)-tData.G'*x0(1:nn)+tData.Vs; % Total energy in structure
    [t,y] = ode45m(@lagTensegrityDynamics,simTime,x0,options,tData);
elseif tData.Correction == 2 % No constraint correction 
    tData.Energy = 0.5*x0(nn+1:end-1)'*tData.M*x0(nn+1:end-1)-tData.G'*x0(1:nn)+tData.Vs; % Total energy in structure
    [t,y] = ode45(@lagTensegrityDynamics,simTime,x0,options,tData);
elseif tData.Correction == 3 % No constraint correction 
    [t,y] = ode15s(@lagTensegrityDynamics_flex,simTime,x0,options,tData);
end   
end