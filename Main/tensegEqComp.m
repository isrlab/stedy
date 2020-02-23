function tData = tensegEqComp(x,tData)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
% 
% Function to find the equilibrium force densities in the strings at t=0, 
% given initial position and minimum force densities in the strings.  
% 
% INPUTS: [x0, x, tData]
% 
% x0: Vector with initial nodal coordinates stacked on top of initial nodal
% velocities. 
% 
% x: Vector with only the initial nodal coordinates. 
% 
% tData.minforce: Scalar field specifiying lower bound for the force densities 
% in the linear programming problem 
%
% NOTE: In case external forces are present during simulation but not
% initially, the user is required to define a vector called 'extF' of
% zeros(size(N)) in another file called extF_eq.m
%
% NEW FIELDS CREATED in tData: [sigmaEq, Energy]
%
% sigmaEq: Vector of output force densities in the strings from the linear
% programming problem. 
%
% Energy: Total Energy in the structure at t=0.
% 

%% Construct gradR
ns = numel(x)/2; % No. of position/velocity variables
q = x(1:ns); % Position vector
qd = x(ns+1:end); % Velocity vector

nLC = size(tData.Lin.A,1); % Number of linear constraints
% nNLC = numel(tData.NLin); % Number of non-linear constraints
nConstr = nLC;% + nNLC;

gradLin = tData.Lin.A; % linear constraint,

% gradNLin = [];
% 
% for i=1:nNLC
%     gradNLin = [gradNLin;x'*tData.NLin(i).Mat];
% end

gradR = gradLin;%;gradNLin];

%% Setting up linear programming problem
% Equilibrium equation for dynamics with compressible bars 
% (eqn. 56) with qdd terms removed.

% Bars
Mdot = zeros(size(tData.M));
Mf = Mdot;
Mq = Mdot;
Xq = zeros(size(tData.M));
beq = zeros(size(Xq.'*q));

bars = tData.bars;
for k=1:tData.nBar
    Xk = tData.listX{k};
    X = tData.X;
    mbk = bars.listM(k); % Mass of bar k
    Ibk = bars.listI(k); % Inertia of bar k
    bk = Xk*q; % Bar k in vector form
    rk = bars.r(k); % Radius of bar k
    lbk = norm(Xk*q); % Length of bar k
    
    lbkdot = bk.'*Xk*qd/lbk;
    rkdot = -bars.nu(k)*rk*lbkdot/lbk;
    Ibkdot = mbk/12*(6*rk*rkdot + 2*lbk*lbkdot);
    
    
    tempMdot = -2*Ibk/lbk^3*(X{k})*lbkdot + Ibkdot/(lbk^2)*(X{k});
    tempMf = Ibk*X{k}*lbkdot/lbk^3;
    tempMq = -(Ibkdot*(X{k})*lbkdot/lbk^3 ...
             - 3*(Ibk*(X{k}))/(lbk^4)*(lbkdot^2) ...
             + (qd.'*(X{k})*q/lbk - q.'*(X{k})*qd*lbkdot/lbk^2)*(Ibk*(X{k}))/lbk^3);
    
    Mdot = Mdot + tempMdot;
    Mf = Mf + tempMf;
    Mq = Mq + tempMq;

end
Mqd = Mdot - Mf;

% Strings with dampers
Fd = zeros(ns,1);
Y = tData.Y;
for k=1:tData.nStr
    lsk = q.'*Y{k}*q;
%     Yk = tData.listY{k};
    c = tData.damper(k);
    lskdot = q.'*Y{k}*qd;

    if isfield(tData,'damper') % If dampers present in structure
        Fd = Fd - c*lskdot*Y{k}*q/lsk;
    end
end

%% bars = tData.bars;
% alpha = cell(1,tData.nBar);
% lb = zeros(tData.nBar,1);
% for k=1:tData.nBar
%     Xk = tData.listX{k};
%     lb(k) = norm(Xk*x); % Length of bar k
%     Ak = pi*bars.r(k)^2;
%     alpha{k} = Ak*bars.E(k)*(Xk.'*Xk)*x;
%     beq = beq + alpha{k}.'*x/lb(k);
% %     lbk0 = bars.L0(k); % Rest length of bar k
% %     Kbk = bars.listK(k); % Stiffness of each bar
% %     
% %     tempXq = Kbk*((Xk.'*Xk) - (Xk.'*Xk)*lbk0/lbk);
% %     
% %     Xq = Xq + tempXq;
% end
% alphaK = cell2mat(alpha);

%% Solving 
Y = cell2mat(tData.Y);
X = cell2mat(tData.X);
L1 = reshape(Y'*q, ns, tData.nStr);
L2 = reshape(X'*q, ns, tData.nBar); % - (minus) for compression
Aeq = [L1 L2 -gradR'];
% beq = Xq.'*x;

if(tData.F)
    run('extF_eq');
    beq = beq + tData.G + extF - Mq.'*q - Mqd.'*qd + Fd;
else beq = beq + tData.G + Fd;
end

% LB = [tData.minforce*ones(tData.nStr,1); tData.minforce*ones(tData.nBar,1);...
%         -Inf*ones(nConstr,1)];%1./lb
% UB = [Inf*ones(tData.nStr + nConstr + tData.nBar,1)];%1./lb];

LB = [tData.minforce*ones(tData.nStr,1); -Inf*ones(tData.nBar,1);...
        -Inf*ones(nConstr,1)];    

UB = [Inf*ones(tData.nStr,1); -tData.minforce*ones(tData.nBar,1);...
         Inf*ones(nConstr,1)];%1./lb];
options = optimoptions('linprog','Display','iter');
[SigLamb,fval,exitflag] = linprog([],[],[],Aeq,beq,LB,UB,[],options);

if(exitflag ~= 1)
    msg = 'Optimization failed. Please check input parameters again.';
    error(msg);
end

sigma = SigLamb(1:tData.nStr);
% rlBar = 1./SigLamb(tData.nStr+1:tData.nStr+tData.nBar);
psi = SigLamb(tData.nStr+1:tData.nStr+tData.nBar);
lamb = SigLamb(tData.nStr+ tData.nBar+1:end);
textSprFD = sprintf('%f \n',sigma); 
textBarFD = sprintf('%f \n',psi); 
fprintf('Force Densities in Strings at Equilibrium: \n %s \n',textSprFD)
fprintf('Force Densities in Bars at Equilibrium: \n %s \n',textBarFD)


lhs = kron(sigma.',eye(ns))*Y.'*q - kron(psi.',eye(ns))*X.'*q ...
        - gradR.'*lamb;
rhs = beq;    
%% Computing Total Energy in Structure at Equilibrium
Vs = 0; % Potential Energy in Strings
if(tData.nStr>0)
    S = tData.N*tData.Cs'; % Strings
    for k=1:tData.nStr        
        s = S(:,k);
        L1 = norm(s); % Initial length between string nodes
        K(k) = tData.strings.E(k)*pi*(tData.strings.r(k))^2/L1 + sigma(k); % Stiffness of string
        Lk(k) = L1*(1-sigma(k)/K(k));
        if Lk(k) < L1
            Vs = Vs+0.5*K(k)*(L1-Lk(k))^2;
        end
    end 
    tData.Lk = Lk;
    tData.K = K;
end 
tData.Vs = Vs;
tData.sigmaEq = sigma; % Force densities in springs at equilibrium
tData.psiEq = psi; % Force densities in bars at equilibrium
end