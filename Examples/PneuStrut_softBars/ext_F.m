% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

%% 3 strut
% extF = 30*[zeros(9,1);sin(t);zeros(4,1);sin(t);0;0;sin(t)]; 
% Sinusoidal external force during simulation

%% Ball
% extF = 30*[zeros(9,1);sin(t);0;0;zeros(6,1);0;sin(t);0;0;0;sin(t);zeros(15,1)]; 
% Sinusoidal external force during simulation

%% Tbar
extF = [zeros(8,1);100*(sin(10*t)+sin(20*t)+cos(15*t)+cos(50*t))/4;zeros(3,1)];