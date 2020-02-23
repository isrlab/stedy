function dxdot = linearlagTensegrityDynamics_flex(t,x,tData)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
% 
% Main ODE function generating the algebraic differential equations to be
% solved using ODE45 and the novel constraint correction method. 
% 

ns = (numel(x)-1)/2; % No. of position/velocity variables
q0 = tData.N(:);
qd0 = zeros(ns,1);
q = x(1:ns) + q0; % Position vector
qd = x(ns+1:end-1) + qd0; % Velocity vector

% q = x(1:ns);
% qd = x(ns+1:end-1);

% Cables
if(tData.nStr>0) % If strings present in structure
    sig_k = zeros(tData.nStr,1); % Force Densities
    for i=1:tData.nStr
        sk = tData.listY{i}*q;
        dsk = tData.listY{i}*qd;
        if tData.Lk(i)<norm(sk)
            sig_k(i) = tData.K(i)*(1 - tData.Lk(i)/norm(sk));
        end
    end
end

% Bars

bars = tData.bars;
psi_k = zeros(tData.nBar,1);
for k=1:tData.nBar
    
    Xk = tData.listX{k};
    
    lbk = norm(Xk*q); % Length of bar k
    lbk0 = bars.L0(k); % Rest length of bar k
    Kbk = bars.listK(k); % Stiffness of each bar
    
    psi_k(k) = Kbk*(1 - lbk0/lbk);
    
end

A = tData.sysSS.A; B = tData.sysSS.B;
dq = q - q0; dqd = qd - qd0;
dsig_k = sig_k - tData.sigmaEq;
dpsi_k = psi_k - tData.psiEq;
dxdot = A*[dq;dqd] + B*[dsig_k;dpsi_k];
% dxdot = A*[q;qd] + B*[sig_k;psi_k];
power = (extF'+Fd')*qd;


dxdot = [dxdot;power];