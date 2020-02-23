function [sysSS,Bf] = linSysCompDescriptor(x,tData)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.
% 
% Linearization with compressible bars (Descriptor form)
% Not using closed-form expression directly. 

ns = (numel(x)-1)/2; % No. of position/velocity variables
q = x(1:ns); % Position vector
qd = x(ns+1:end-1); % Velocity vector

Forces = tData.G; % Initializing with gravitational force
Fd = zeros(ns,1); % Force due to dampers
dFd_dq = zeros(ns,ns);
dFd_dqd = zeros(ns,ns);
Y = tData.Y;
% Cables
if(tData.nStr>0) % If strings present in structure
    for k=1:tData.nStr
        lsk = q.'*Y{k}*q;
        
        c = tData.damper(k);
        lskdot = q.'*Y{k}*qd;
        dlsk_dq = q.'*Y{k}/lsk ;
        dlskd_dq = qd.'*Y{k}/lsk - q.'*Y{k}*qd/lsk^2*dlsk_dq;
        dlskd_dqd = q.'*Y{k}/lsk;

        if isfield(tData,'damper') % If dampers present in structure
            dFd_dq = -c*(Y{k}*q/lsk*dlskd_dq + lskdot*Y{k}/lsk ...
                        - lskdot*Y{k}*q/lsk^2*dlsk_dq) ...
                     + dFd_dq;
            dFd_dqd = -c*Y{k}*q/lsk*dlskd_dqd + dFd_dqd;
        end
    end
end

sig_k = tData.sigmaEq; % FD in strings at Eq
psi_k = tData.psiEq; % FD in bars at Eq
% Bars
% Mqdd, Mqd, Mq
% Mqd = Mdot - Mf
Mqdd = tData.M;
Mdot = zeros(size(tData.M));
Mf = Mdot;
Mq = Mdot;
% dMqdT_dqd = zeros(ns^2,ns);
dMqd_dqdxqd = zeros(size(tData.M));
dMqd_dqxqd = zeros(size(tData.M));
dMq_dqdxq = zeros(size(tData.M));
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
    
    dlbk_dq = q.'*X{k}/lbk;
    dlbkd_dq = qd.'*X{k}/lbk - qd.'*X{k}*q/lbk^2*dlbk_dq;
    dlbkd_dqd = dlbk_dq;
    dIbk_dq = mbk*lbk/6*dlbk_dq;
    drkd_dq = -bars.nu(k)*rk/lbk*dlbkd_dq + ...
                bars.nu(k)*rk*lbkdot/lbk^2*dlbk_dq;
    d_dIbkdt_dq = mbk/12*(6*rk*drkd_dq + 2*lbkdot*dlbk_dq...
                        + 2*lbk*dlbkd_dq);
    d_dIbkdt_dqd = mbk/12*(-6*rk^2*bars.nu(k)/lbk^2 + 2)*dlbk_dq;

    tempMdot = -2*Ibk/lbk^3*(X{k})*lbkdot + Ibkdot/(lbk^2)*(X{k});
    tempMf = Ibk*X{k}*lbkdot/lbk^3;
    tempMq = -(Ibkdot*(X{k})*lbkdot/lbk^3 ...
             - 3*(Ibk*(X{k}))/(lbk^4)*(lbkdot^2) ...
             + (qd.'*(X{k})*q/lbk - q.'*(X{k})*qd*lbkdot/lbk^2)*(Ibk*(X{k}))/lbk^3);
    tempMqdd = -Ibk*(X{k})*(q*q.')*(X{k})/lbk^4;
    tempXq = Kbk*((X{k}) - (X{k})*lbk0/lbk);
    
%     tempdMqd_dqd = vec(- 3*Ibk/lbk^3*(X{k}) + d_dIbkdt_dqd/lbk^2*(X{k}));
    tempdMqd_dqdxqd = -3*Ibk/lbk^3*X{k}*qd*dlbkd_dqd ...
                      + X{k}/lbk^2*qd*d_dIbkdt_dqd;
    
    tempdMqd_dqxqd = -3*Ibk*X{k}*qd/lbk^3*dlbkd_dq - ...
                    - 3*Ibk*lbkdot*X{k}*qd/lbk^4*dlbk_dq...
                    + X{k}*qd/lbk^2*d_dIbkdt_dq...
                    -2*Ibkdot*X{k}*qd/lbk^3*dlbk_dq...
                    -3*lbkdot*X{k}*qd/lbk^3*dIbk_dq;
                
    tempdMq_dqdxq = X{k}*q*lbkdot/lbk^3*d_dIbkdt_dqd ...
                    + Ibkdot*X{k}*q/lbk^3*dlbkd_dqd ...
                    - 6*Ibk*X{k}*q*lbkdot/lbk^4*dlbkd_dqd ...
                    + (Ibk*X{k}/lbk^3)*(2*X{k}*q*qd.'/lbk ...
                    - q*q.'*X{k}*lbkdot/lbk^2 ...
                    - q*q.'*X{k}*qd/lbk^2*dlbkd_dqd);
    
    Mdot = Mdot + tempMdot;
    Mf = Mf + tempMf;
    Mq = Mq + tempMq;
    Mqdd = Mqdd + tempMqdd;
    Xq = Xq + tempXq;
    dMqd_dqdxqd = dMqd_dqdxqd + tempdMqd_dqdxqd;
    dMqd_dqxqd = dMqd_dqxqd + tempdMqd_dqxqd;
    dMq_dqdxq = dMq_dqdxq + tempdMq_dqdxq;
end
Mqd = Mdot - Mf;

% Construct various matrices
nLC = size(tData.Lin.A,1); % Number of linear constraints
% nConstr = nLC;
gradLin = tData.Lin.A; % Jacobian of the linear constraints,

Rq = gradLin;

YTq = cell2mat(tData.Y).'*q;
Yhat = reshape(YTq,[ns,tData.nStr]);
XTq = cell2mat(tData.X).'*q;
Xhat = reshape(XTq,[ns,tData.nBar]);

Ksys = genKsys(x,tData);
% xi3 = Cab_En - Mq.'*q - Xq.'*q - Mqd.'*qd + tData.G + extF + Fd;
dxi3_dq = -(Mq.' + kron(tData.sigmaEq.',eye(ns))*cell2mat(Y).' ...
         + kron(tData.psiEq.',eye(ns))*cell2mat(X).') ...
         - dMqd_dqxqd + dFd_dq; 
% dxi3_dq = -(Mq.' + Ksys) - dMqd_dqxqd + dFd_dq; 
dxi3_dqd = -dMq_dqdxq - Mqd.' - dMqd_dqdxqd + dFd_dqd;
dxi3_dsigma = -Yhat;
dxi3_dpsi = -Xhat;
dxi3_df = eye(ns,ns);

dxi3_dq0 = -Mq.' - dMqd_dqxqd + dFd_dq;

M1 = [Mqdd -Rq.'; -Rq zeros(nLC,nLC)];
M_alpha = inv(M1);
M_beta = M_alpha(1:ns,1:ns);

% EB = eye(ns) - Rq.'*pinv(Rq.');
% FC = eye(ns) - pinv(Rq)*Rq;
% r1 = rank(Rq.')
% r2 = rank(Rq)
% r3 = rank(EB*Mqdd*FC)
% szM1 = size(M1)
% rM1 = rank(M1)

A = [zeros(ns,ns) eye(ns,ns); 
            M_beta*dxi3_dq M_beta*dxi3_dqd];
        
A0 = [zeros(ns,ns) eye(ns,ns);
      M_beta*dxi3_dq0 M_beta*dxi3_dqd];
eA0 = eig(A0);
% scatter(real(eA0),imag(eA0));
Bu = [zeros(ns,tData.nStr) zeros(ns,tData.nBar); 
             M_beta*dxi3_dsigma M_beta*dxi3_dpsi];
         
Bf = [zeros(ns,ns); M_beta*dxi3_df];

C = eye(size(A));
D = zeros(size(Bu));
sysSS = ss(A, Bu, C, D);
end