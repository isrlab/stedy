function xdot = lagTensegrityDynamics_flex(t,x,tData)
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
psi_k = zeros(tData.nBar,1);
for k=1:tData.nBar
    X = tData.X;
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
    psi_k(k) = Kbk*(1 - lbk0/lbk);
    tempXq = Kbk*(X{k} - X{k}*lbk0/lbk);
    
    Mdot = Mdot + tempMdot;
    Mf = Mf + tempMf;
    Mq = Mq + tempMq;
    Mqdd = Mqdd + tempMqdd;
    Xq = Xq + tempXq;
end
Bar_En = -kron(psi_k',eye(ns))*(cell2mat(tData.X))'*q;
Mqd = Mdot - Mf;


% % Bars Potential Energy
%     sig_k_bar = zeros(tData.nBar,1); % Force Densities
%     for i=1:tData.nBar
%         sk = tData.listX{i}*q;
%         dsk = tData.listX{i}*qd;
%         sig_k_bar(i) = tData.bars.K(i)*(1 - tData.bars.L0(i)/norm(sk));
%     end
%     Cab_En = Cab_En-kron(sig_k_bar',eye(ns))*(cell2mat(tData.NLin))'*q;

% Construct various matrices
nLC = size(tData.Lin.A,1); % Number of linear constraints
nConstr = nLC;
gradLin = tData.Lin.A; % Jacobian of the linear constraints,
hessLin = zeros(nLC,1); % Hessian of the linear constraints,

gradR = [gradLin;];
hessR = [hessLin;];

M = [Mqdd -gradR';
     -gradR zeros(nConstr,nConstr)];
    
% Check for external force
if tData.F == 1
    run('ext_F');
else
    extF = zeros(ns,1);
end
    
if(tData.nStr>0) % If strings present in structure
    F1_qdd = Cab_En - Mq.'*q + Bar_En - Mqd.'*qd + tData.G + extF + Fd;
    F1 = [F1_qdd;hessR]; % No gravity
else
    F1 = [tData.G+extF + Fd;hessR];
end


xSol = M\F1; % Solve the linear equation
qdd = xSol(1:ns);
power = (extF'+Fd')*qd;


xdot = [qd;qdd;power];