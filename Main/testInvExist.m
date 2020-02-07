function xdot = lagTensegrityDynamics_flex(x,tData)
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

% Bars
% Mqdd, Mqd, Mq
% Mqd = Mdot - Mf
Mqdd = tData.M;
Mdot = zeros(size(tData.M));
Mf = Mdot;
Mq = Mdot;
Xq = Mdot;
bars = tData.bars;
for k=1:tData.nBar
    Xk = tData.listX{k};
    mbk = bars.listM(k); % Mass of bar k
    Ibk = bars.listI(k); % Inertia of bar k
    bk = Xk*q; % Bar k in vector form
    rk = bars.r(k); % Radius of bar k
    lbk = norm(Xk*q); % Length of bar k
    lbk0 = bars.L0(k); % Rest length of bar k
    Kbk = bars.listK(k); % Stiffness of each bar
    
    lbkdot = bk.'*Xk*qd/lbk;
    rkdot = -bars.nu(k)*rk*lbkdot/lbk;
    Ibkdot = mbk/12*(6*rk*rkdot + 2*lbk*lbkdot);
    
    tempMdot = -2*Ibk/lbk^3*(Xk.'*Xk)*lbkdot + Ibkdot/(lbk^2)*(Xk.'*Xk);
    tempMf = Ibk*Xk.'*Xk*lbkdot/lbk^3;
    tempMq = -(Ibkdot*(Xk.'*Xk)*lbkdot/lbk^3 ...
             - 3*(Ibk*(Xk.'*Xk))/(lbk^4)*(lbkdot^2) ...
             + (qd.'*(Xk.'*Xk)*q/lbk - q.'*(Xk.'*Xk)*qd*lbkdot/lbk^2)*(Ibk*(Xk.'*Xk))/lbk^3);
    tempMqdd = -Ibk*(Xk.'*Xk)*q*q.'*(Xk.'*Xk)/lbk^4;
    tempXq = Kbk*((Xk.'*Xk) - (Xk.'*Xk)*lbk0/lbk);
    
    Mdot = Mdot + tempMdot;
    Mf = Mf + tempMf;
    Mq = Mq + tempMq;
    Mqdd = Mqdd + tempMqdd;
    Xq = Xq + tempXq;
end
Mqd = Mdot - Mf;

% Construct various matrices
nLC = size(tData.Lin.A,1); % Number of linear constraints
nConstr = nLC;
gradLin = tData.Lin.A; % Jacobian of the linear constraints,

Rq = [gradLin;];

M = [Mqdd -Rq';
     -Rq zeros(nConstr,nConstr)];
sz = size(M)
r = rank(M)
szMqdd = size(Mqdd)
rMqdd = rank(Mqdd)
Minv = inv(M);    
% Construct C matrix (inv of M)
C1 = null(Rq);
C2 = (-Rq.')*inv(Rq*Rq.');
C3 = inv(Rq*Rq.')*Rq*(Mqdd*C1 - eye(size(Mqdd*C1)));
C4 = inv(Rq*Rq.')*Rq*(Mqdd*C2);
C = [C1 C2; C3 C4];

C1_new = [C1 2*C1 3*C1];
C3_new = inv(Rq*Rq.')*Rq*...
            (Mqdd*C1_new - eye(size(Mqdd*C1_new)));
C_new = [C1_new C2; C3_new C4];

C1true = Minv(1:ns,1:ns);
C2true = Minv(1:ns,ns+1:end);
C3true = Minv(ns+1:end,1:ns);
C4true = Minv(ns+1:end,ns+1:end);
c = 0;