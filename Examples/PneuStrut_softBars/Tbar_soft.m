% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
% 
% To compare Tbar w/o correction against Tbar_flex.
clc; clear;
close all;
%% Define initial nodal coordinates and connectivity matrices
theta = pi/4; % To help describe angular position of any node

N = [0  5*cos(theta)  5*cos(theta)              0;%    5*sin(theta);
     0            0              0              0;%    0;
     0            0   5*sin(theta)   5*sin(theta);];%   10*sin(theta);];
 
% N = [0  6*cos(theta)  6*cos(theta)     cos(theta);
%      0            0              0              0;
%      0            0   6*sin(theta)   5*sin(theta);]; 

% fixedNodes = zeros(size(N));
% fixedNodes = [1 1 0; 1 1 0; 0 1 0; 0 1 0;]';% 0 1 0]'; % To identify which coordinates of a node are fixed: 1-fixed, 0-unfixed
                                            % Same size as N

fixedNodes = [1 1 1; 1 1 1; 0 1 0; 0 1 0;]'; % 

% Connectivity Matrices
% Cb an Cs have been defined independently. 
Cb = [-1  0  1  0;%  0; 
       0 -1  0  1;];%  0;
%        0 0  0 -1 1;];
       
Cs = [-1  1  0  0;% 0;
       0 -1  1  0;% 0;
       0  0 -1  1;% 0;
       1  0  0 -1;];% 0;]; 

tData = tensegConstruct(N,Cb,Cs,fixedNodes);

% Vector of Initial Nodal Positions and Velocities
% v0 = [0 0 0 0;
%       0 0 0 0;
%       1 1 1 1;];
% x0 = [N(:); v0(:)];  
x0 = [N(:); 0*N(:)];

%% Material and Simulation Environment properties
% Strings
strings.r = 1/1000*ones(1,tData.nStr); % Radius of strings
strings.E = 7e7*ones(1,tData.nStr); % Young's Modulus of Nylon
strings.rLP = ones(1,tData.nStr);  % Rest lengths of the strings: 0.7 means 70% 
% strings.rLP = [1 0.7 1 0.7];  % Rest lengths of the strings: 0.7 means 70% 

% Bars
barsSoft.r = 0.05*ones(1,tData.nBar); % Radius of bars
barsSoft.rho = 960*ones(1,tData.nBar); % Density of bars
barsSoft.nu = 0.46*ones(1,tData.nBar); % Poisson's ratio of bars (HDPE)
barsSoft.E = 1e9*ones(1,tData.nBar); % Young's modulus of bars (HDPE)

barsMetal.r = 0.05*ones(1,tData.nBar); % Radius of bars
barsMetal.rho = 2700*ones(1,tData.nBar); % Density of bars
barsMetal.nu = 0.30*ones(1,tData.nBar); % Poisson's ratio of bars (aluminium)
barsMetal.E = 200e9*ones(1,tData.nBar); % Young's modulus of bars

% Point Masses
Mp = ones(1,tData.nPm); % All point masses initialised with a mass of 1

g = [0;0;0]; % Gravity

tData = tensegGenMat(tData,barsSoft,strings,Mp,g);
tDataRigid = tensegGenMat(tData,barsMetal,strings,Mp,g);

% The structure is in equilibrium in the initial position. Hence, there
% is no need to call the function tensegEq.

%% Simulation Inputs
tData.F = 0; % If 1, external forces are present in structure, not if 0.
tDataRigid.F = 0;
tData.damper = ones(1,tData.nStr); % All strings initialised with dampers whose damping coefficient is 1. 
tDataRigid.damper = ones(1,tData.nStr); % All strings initialised with dampers whose damping coefficient is 1. 
tData.minforce = 0; % Lower bound for force densities in the strings
tDataRigid.minforce = 0; % Lower bound for force densities in the strings
tData = tensegEqComp(x0,tData);
tData.Compressible = 1;
tDataRigid = tensegEqComp(x0,tDataRigid);
tDataRigid.Compressible = 1;
%% Checking genOptK
ns = numel(N);
nC = nnz(fixedNodes);
Y = cell2mat(tData.Y);
X = cell2mat(tData.X);
% Ksys = kron(tData.sigmaEq.',eye(ns))*Y.' ...
%         + kron(tData.psiEq.',eye(ns))*X.';

Ksys = genKsys(x0,tData);
% KsysRigid = kron(tDataRigid.sigmaEq.',eye(ns))*Y.'; 
KsysRigid = genKsys(x0,tDataRigid);    
%% Removing constrained rows and columns
% 1st, 2nd nodes are constrained in x,y,z.
% 3rd and 4th nodes are constrained in y. 
% Remaining rows/cols:7,9,10,12/
W = zeros(ns,ns-nC);
W(7,1) = 1; W(9,2) = 1; W(10,3) = 1; W(12,4) = 1;
Kred = W.'*Ksys*W;
KredRigid = W.'*KsysRigid*W;
% % % % % % % % Krem = Ksys;
% % % % % % % % Krem(1,:) = []; Krem(:,1) = [];
% % % % % % % % Krem(1,:) = []; Krem(:,1) = [];
% % % % % % % % Krem(1,:) = []; Krem(:,1) = [];
% % % % % % % % Krem(1,:) = []; Krem(:,1) = [];
% % % % % % % % Krem(1,:) = []; Krem(:,1) = [];
% % % % % % % % Krem(1,:) = []; Krem(:,1) = [];
% % % % % % % % Krem(2,:) = []; Krem(:,2) = [];
% % % % % % % % Krem(4,:) = []; Krem(:,4) = [];

%% Deflection
s = svd(Kred);    
[U,~,V] = svd(Kred);
sRigid = svd(KredRigid);
[Urigid,~,Vrigid] = svd(KredRigid);

% ind = find(s > 1e-12);
OneByK = zeros(ns-nC,ns-nC);
OneByKRigid = zeros(ns-nC,ns-nC);
for i = 1:length(s)
   OneByK = OneByK + 1/s(i)*V(:,i)*(U(:,i).');
   OneByKRigid = OneByKRigid + ...
                   1/sRigid(i)*Vrigid(:,i)*(Urigid(:,i).');
end
Fext = ones(ns-nC,1);
Fext(1) = 0; Fext(3) = 0; Fext(4) = 0;
% Fext(3) = 0; Fext(6) = 0; Fext(9) = 10; Fext(12) = 10;
% Fext(3) = 0; Fext(6) = 0; Fext(8) = 10; Fext(11) = 10;
% Fext(3) = 0; Fext(6) = 0; Fext(7) = 10; Fext(10) = 10;
deltx = OneByK*Fext
deltxRigid = OneByKRigid*Fext

% KsysNew = Ksys; 
% % KsysNew(3,:) = 2*Ksys(3,:); KsysNew(:,3) = 2*Ksys(:,3); 
% % KsysNew(6,:) = 2*Ksys(6,:); KsysNew(:,6) = 2*Ksys(:,6);
% KsysNew(9,:) = 2*Ksys(9,:); KsysNew(:,9) = 2*Ksys(:,9);
% KsysNew(12,:) = 2*Ksys(12,:); KsysNew(:,12) = 2*Ksys(:,12);
% % KsysNew(10:12,10:12) = Ksys(4:6,4:6);
% % KsysNew = Ksys + (1e-6)*rand(size(Ksys));
% % KsysNew = 10*Ksys;
% [sigOpt, psiOpt] = genOptK(tData,KsysNew);
% Kopt = kron(sigOpt.',eye(ns))*Y.' ...
%         + kron(psiOpt.',eye(ns))*X.';
% figure(); imagesc(abs(Kopt - KsysNew)); colorbar;    
% err = norm(sigOpt - tData.sigmaEq) + norm(psiOpt - tData.psiEq);
%% Final Simulation
tEnd = 10; % Simulation End Time

x0 = [x0;0]; % Initial Condition - [Position; Velocity; Energy];

options = odeset('RelTol',1e-10,'AbsTol',1e-10,'Refine',1);

% [simTime,tInt] = tensegSimTime(options,tEnd);

tData.Correction = 3; % Compressible Bar 
% [tFlex,yFlex] = tensegSim(x0,simTime,tData,options);

[sysSS, Bf] = linSysCompDescriptor(x0, tData);
tData.sysSS = sysSS;

%% Minreal
minSys = minreal(sysSS);%,10);
Amin = minSys.A;
szAmin = size(Amin)
rAmin = rank(Amin)
% 
% % [balSys,g] = balreal(sysSS)
% % rSys = balred(sysSS,szAmin(1))
% bodeplot(sysSS,minSys,'r--')
% Co = ctrb(minSys);
% unco = length(minSys.A) - rank(Co)

% Using ctrbf to get the controllable subspace
% [Abar,Bbar,Cbar,T,k] = ctrbf(sysSS.A,sysSS.B,sysSS.C);
% nCntrlSt = sum(k);
% Ac = Abar(end - nCntrlSt + 1: end, end - nCntrlSt + 1: end);
% Bc = Bbar(end - nCntrlSt + 1: end, :);
% nsCntrl = size(Ac,1);
% nuCntrl = size(Bc,2);
% ns = numel(N);
% nu = tData.nStr + tData.nBar;
% nsMin = size(minSys.A,1); nuMin = size(minSys.B,2);

%% H2
% % % % Bw = T*Bf; 
% % % % Bw= Bw(end - nCntrlSt + 1:end,:);
% % % % Cz = T(end - nCntrlSt + 1:end,:).';
% % % % Du = zeros(2*ns,nuCntrl);
% % % % gam = 1;
% % % % cvx_begin sdp 
% % % %     variable X(nsCntrl, nsCntrl) symmetric
% % % %     variable W(2*ns, 2*ns) symmetric
% % % %     variable Z(nuCntrl, nsCntrl)
% % % %     
% % % %     [Ac Bc]*[X; Z]  + [X Z']*[Ac'; Bc'] + Bw*Bw' <= 0
% % % %     [W Cz*X; X*Cz' X] >= 0
% % % %     trace(W) < gam
% % % %     
% % % %     minimize trace(W)
% % % % cvx_end
% % % % h2K = Z*inv(X);
% % % % KtrH2 = [zeros(nuCntrl, length(Abar) - nsCntrl) h2K]*T;
% % % % sysH2 =  ss(sysSS.A - sysSS.B*KtrH2, Bf,...
% % % %                     sysSS.C, zeros(2*ns,ns));   
    
%% LQR
% Q = 1e-2*eye(nsCntrl);
% R = eye(nuCntrl);
% 
% [Klqr,S,E] = lqr(Ac,Bc,Q,R);
% % % % % [KlqrDist,S,E] = lqr(Ac,Bc,Cz.'*Cz,Du.'*Du);
% Kaug = [zeros(nuCntrl, length(Abar) - nsCntrl) Klqr];
% Ktrlqr = Kaug*T;
% syslqr = ss(sysSS.A - sysSS.B*Ktrlqr,zeros(size(sysSS.B)),...
%                     sysSS.C, sysSS.D);
% % % % % sqslqrDist = ss(sysSS.A - sysSS.B*Ktrlqr,Bf,...
% % % %                     sysSS.C, sysSS.D);
%% Sim                
% tSim = (0:0.01:100)';
% x0Sim = zeros(2*ns,1);
% x0Sim(9) = 0.1*rand(1);
% x0Sim(12) = 0.05*rand(1);
% % uSim = 0.1*rand(12,length(tSim));
% [ySim, tSim, xSim] =  ...
%     lsim(syslqr,zeros(nuCntrl,length(tSim)),tSim,x0Sim);
% % [ySim, tSim, xSim] =  ...
% %     lsim(sysH2,uSim,tSim,x0Sim);
% uSim = zeros(length(tSim),nuCntrl);
% for i=1:length(tSim)
%     uSim(i,:) = Ktrlqr*xSim(i,:)';
% end
%% Plot
% close all;
% figure(1); clf;
% subplot(4,1,1); plot(tSim, ySim(:,7));title('Position'); ylabel('3x');
% subplot(4,1,2); plot(tSim, ySim(:,9)); ylabel('3z');
% subplot(4,1,3); plot(tSim, ySim(:,10)); ylabel('4x');
% subplot(4,1,4); plot(tSim, ySim(:,12)); ylabel('4z'); xlabel('Time(s)');
% 
% figure(2); clf;
% subplot(4,1,1); plot(tSim, ySim(:,19)); title('Velocity'); ylabel('3x');
% subplot(4,1,2); plot(tSim, ySim(:,21)); ylabel('3z');
% subplot(4,1,3); plot(tSim, ySim(:,22)); ylabel('4x');
% subplot(4,1,4); plot(tSim, ySim(:,24)); ylabel('4z'); xlabel('Time(s)');
% 
% figure(3); clf;
% subplot(4,1,1); plot(tSim, uSim(:,1));title('Cable FD'); ylabel('s1');
% subplot(4,1,2); plot(tSim, uSim(:,2)); ylabel('s2');
% subplot(4,1,3); plot(tSim, uSim(:,3)); ylabel('s3');
% subplot(4,1,4); plot(tSim, uSim(:,4)); ylabel('s4'); xlabel('Time(s)');
% 
% figure(4); clf;
% subplot(2,1,1); plot(tSim, uSim(:,5));title('Bar FD'); ylabel('b1');
% subplot(2,1,2); plot(tSim, uSim(:,6)); ylabel('b2');xlabel('Time(s)');

% figure();clf
% plot(tSim, ySim(:,1));
% [K,S,E] = lqr(minSys,Q,R,Nlqr);
%% Checking Linear Model against Nonlinear
x0delt = zeros(size(x0));%(1:end-1);
T = 0:0.01:0.2;
x0delt(9) = 0.1*rand(1);
[tSim,ySim] = ode15s(@linearlagTensegrityDynamics_flex,T,x0delt,options,tData);
% % [ysim,tsim,xsim] = lsim(sysSS,zeros(length(T),6),T,x0delt);
% % plot(tSim,ySim(:,12));
[tFlex,yFlex] = tensegSim(x0delt+x0,T,tData,options);
for i=1:length(tSim)
    ySim(i,:) = ySim(i,:) + x0(1:end)';
end
%% Plotting 
% Plot Output Trajectories
plotMotion(tFlex, yFlex, tData);
plotMotion(tSim, ySim, tData);

% Overlay Motion Plots
plotCompTbar_flex(tSim,ySim,tFlex,yFlex,tData);

% print(figure(3),'compTBarFlex_MotionNode3','-depsc');
% print(figure(4),'compTBarFlex_MotionNode4','-depsc');

