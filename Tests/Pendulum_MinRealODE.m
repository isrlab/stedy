function qDot = Pendulum_MinRealODE(t,q,tData)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
% 
% ODE file for the minimum realization simulation of a single pendulum,
% conventionally accepted as the ground truth in dynamics circles. 
% 
q1 = q(1);
q2 = q(2);

g = 9.806;

A = ...
[ 1,     0;
  0, -25/3];
 
b = ... 
[              q2;
 (5*g*sin(q1))/2];

if(rank(A)<2)
    keyboard
end
qDot = [A\b];
end
