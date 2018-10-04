function xdot = lagTensegrityDynamics(t,x,tData)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
% 
% Main ODE function generating the algebraic differential equations to be
% solved using ODE45 and the novel constraint correction method. 
% 

ns = (numel(x)-1)/2; % No. of position/velocity variables
q = x(1:ns); % Position vector
qd = x(ns+1:end-1); % Velocity vector

Forces = tData.G; % Initializing with gravitational force
Fd = zeros(ns,1);

% Cables
if(tData.nStr>0) % If strings present in structure
    sig_k = zeros(tData.nStr,1); % Force Densities
    for i=1:tData.nStr
        sk = tData.listY{i}*q;
        dsk = tData.listY{i}*qd;
        if tData.Lk(i)<norm(sk)
            sig_k(i) = tData.K(i)*(1 - tData.Lk(i)/norm(sk));
            if isfield(tData,'damper') % If dampers present in structure
                Fd = Fd+((-tData.damper(i)*(dsk'*sk)*sk/(sk'*sk))'*tData.listY{i})';
            end
        end
        Forces = Forces - sig_k(i)*tData.Y{i}'*tData.Y{i}*q;
    end
    Cab_En = -kron(sig_k',eye(ns))*(cell2mat(tData.Y))'*q;

end

% Construct various matrices
nLC = size(tData.Lin.A,1); % Number of linear constraints
nNLC = numel(tData.NLin); % Number of non-linear constraints
nConstr = nLC + nNLC;

gradLin = tData.Lin.A; % Jacobian of the linear constraints,
hessLin = zeros(nLC,1); % Hessian of the linear constraints,

gradNLin = []; % Jacobian of the non-linear constraints
hessNLin = []; % Hessian of the non-linear constraints

for i=1:nNLC
    gradNLin = [gradNLin;2*q'*tData.NLin(i).Mat];
    hessNLin = [hessNLin;2*qd'*tData.NLin(i).Mat*qd];

end
gradR = [gradLin;gradNLin];
hessR = [hessLin;hessNLin];

M = [tData.M -gradR';
     -gradR zeros(nConstr,nConstr)];

% Check for external force
if tData.F == 1
    run('ext_F');
else
    extF = zeros(ns,1);
end
    
if(tData.nStr>0) % If strings present in structure
    F1 = [tData.Cab_En + tData.G + extF + Fd;hessR]; % No gravity
else
    F1 = [tData.G+extF + Fd;hessR];
end

xSol = M\F1; % Solve the linear equation
qdd = xSol(1:ns);
power = (extF'+Fd')*qd;

xdot = [qd;qdd;power];