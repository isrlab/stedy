function linSys = linSysComp(x,tData)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
% 
% Linearization with compressible bars.

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
dMqT_dqd = zeros(ns^2,ns);
Xq = Mdot;
bars = tData.bars;

dqqT_dqi = cell(ns,1);
dMqdd_dq = cell(1,ns);
for i=1:ns
    dqqT_dqi{i} = zeros(ns,ns);
    dqqT_dqi{i}(i,:) = q.';
    dqqT_dqi{i}(:,i) = q;
    dqqT_dqi{i}(i,i) = 2*q(i);
end
for i=1:ns
    tempMqdd_dqi = zeros(ns,ns);
    for k=1:tData.nBar
        Xk = tData.listX{k};
        X = tData.X;
        mbk = bars.listM(k); % Mass of bar k
        Ibk = bars.listI(k); % Inertia of bar k
        lbk = norm(Xk*q); % Length of bar k
        
        tempMqdd_dqi = tempMqdd_dqi ...
                       + Ibk*(X{k})/(lbk^4)*dqqT_dqi{i}*(X{k}) ...
                       + mbk/12/(lbk^5)*q(i)*(X{k})*(X{k})*(q*q.')*(X{k}) ...
                       - 2*q(i)*(X{k})*Ibk*(X{k})/(lbk^6)*(q*q.')*(X{k});
        
    end
    dMqdd_dq{i} = tempMqdd_dqi;
end

for k=1:tData.nBar
    Xk = tData.listX{k};
    X = tData.X;
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
    Ibkdd = mbk/12*(-6*rk^2*bars.nu(k)/lbk^2 + 2);
    
    tempMdot = -2*Ibk/lbk^3*(X{k})*lbkdot + Ibkdot/(lbk^2)*(X{k});
    tempMf = Ibk*X{k}*lbkdot/lbk^3;
    tempMq = -(Ibkdot*(X{k})*lbkdot/lbk^3 ...
             - 3*(Ibk*(X{k}))/(lbk^4)*(lbkdot^2) ...
             + (qd.'*(X{k})*q/lbk - q.'*(X{k})*qd*lbkdot/lbk^2)*(Ibk*(X{k}))/lbk^3);
    tempMqdd = -Ibk*(X{k})*(q*q.')*(X{k})/lbk^4;
    tempXq = Kbk*((X{k}) - (X{k})*lbk0/lbk);
    tempdMq_dqd = vec(- Ibk/lbk^4*(X{k}) + Ibkdd/lbk^2*(X{k}));
    
    Mdot = Mdot + tempMdot;
    Mf = Mf + tempMf;
    Mq = Mq + tempMq;
    Mqdd = Mqdd + tempMqdd;
    Xq = Xq + tempXq;
    dMqT_dqd = dMqT_dqd + tempdMq_dqd*q.'*(X{k});
end
Mqd = Mdot - Mf;
dMqT_dqdxq = reshape(dMqT_dqd*q, [ns,ns]);
% Construct various matrices
nLC = size(tData.Lin.A,1); % Number of linear constraints
nConstr = nLC;
gradLin = tData.Lin.A; % Jacobian of the linear constraints,
hessLin = zeros(nLC,1); % Hessian of the linear constraints,

Rq = [gradLin;];
hessR = [hessLin;];
   
% Check for external force
if tData.F == 1
    run('ext_F');
else
    extF = zeros(ns,1);
end

xi3 = Cab_En - Mq.'*q - Xq.'*q - Mqd.'*qd + tData.G + extF + Fd;
dxi3_dq = -(Mq.' + Xq.'+ ...
            kron(sig_k',eye(ns))*(cell2mat(tData.Y))'); 
Mqdd_inv = inv(Mqdd);
Mqdd_tilde = inv(Rq*Mqdd_inv*Rq.');
xi2 = hessR;
xi_tilde = xi2 + Rq*Mqdd_inv*xi3;

dMqddInv_dq = cell(1,ns);
dMqdd_tilde_dq = cell(1,ns);
dxi_tilde_dq = cell(1,ns);
dzeta_dq1 = zeros(ns,ns); % first term in dzeta_dq
for i=1:ns
    dMqddInv_dq{i} = -Mqdd_inv*dMqdd_dq{i}*Mqdd_inv;
    dzeta_dq1 = dzeta_dq1 + ...
                dMqddInv_dq{i}*(xi3 - Rq.'*Mqdd_tilde*xi_tilde);
    dMqdd_tilde_dq{i} = vec(-inv(Mqdd_tilde)*Rq*...
                        dMqddInv_dq{i}*(Rq.')*inv(Mqdd_tilde));
    dxi_tilde_dq{i} = Rq*dMqddInv_dq{i}*xi3;%
end
dMqddInv_dq = cell2mat(dMqddInv_dq);
dMqdd_tilde_dq = cell2mat(dMqdd_tilde_dq);
dxi_tilde_dq = cell2mat(dxi_tilde_dq);
dxi_tilde_dq = dxi_tilde_dq + Rq*Mqdd_inv*dxi3_dq;

YTq = cell2mat(tData.Y).'*q;
Yhat = reshape(YTq,[ns,tData.nStr]);
XTq = cell2mat(tData.X).'*q;
Xhat = reshape(XTq,[ns,tData.nBar]);

dzeta_dq = dzeta_dq1 ...
           + Mqdd_inv*(dxi3_dq - kron(xi_tilde.',Rq.')*dMqdd_tilde_dq ...
           - Rq.'*Mqdd_tilde*dxi_tilde_dq);
dzeta_dqd = Mqdd_inv*(-dMqT_dqdxq - Mqd.' - ...
            Rq.'*Mqdd_tilde*Rq*Mqdd_inv*(-dMqT_dqdxq - Mqd.'));
        
dzeta_dsigma = Mqdd_inv*(-Yhat - Rq.'*Mqdd_tilde*Rq*Mqdd_inv*(-Yhat));
dzeta_dpsi = Mqdd_inv*(-Xhat - Rq.'*Mqdd_tilde*Rq*Mqdd_inv*(-Xhat));
dzeta_df = Mqdd_inv*(eye(ns,ns) - Rq.'*Mqdd_tilde*Rq*Mqdd_inv);

linSys.A = [zeros(ns,ns) eye(ns,ns); dzeta_dq dzeta_dqd];
linSys.Bu = [zeros(ns,tData.nStr) zeros(ns,tData.nBar); dzeta_dsigma dzeta_dpsi];
linSys.Bf = [zeros(ns,ns); dzeta_df];
end