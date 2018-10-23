function [simTime, tInt] = tensegSimTime(options,tEnd)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.

% Function to generate simulation time and output interval time-step for
% simulation and animation. 
% 
% INPUTS: 
% options: ODE45 options. The user is encouraged to look at MATLAB's ode45
% documentation to understand how to provide this input.
% 
% tEnd: Simulation end time
% 
% OUTPUTS:
% tInt: Output interval time-step
% 
% simTime: Simulation time, to be provided as input to tensegSim. If Refine
% option is turned on, simTime is just [0 tSpan], otherwise the function
% prompts the user to enter output interval time-step and simTime shall form a
% vector [0:tInt:tEnd]
% 
if (options.Refine)
    simTime = [0 tEnd];
    tInt = 0.01;
else
    prompt = 'Enter output time-step (e.g. 0.01): ';
    tInt = input(prompt);
    simTime = 0:tInt:tEnd;
end

end