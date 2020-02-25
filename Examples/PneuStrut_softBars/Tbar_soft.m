% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
% 
% To compare Tbar w/o correction against Tbar_flex.
clc; clear;
% close all;
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
strings.E = 2E9*ones(1,tData.nStr); % Young's Modulus of Nylon
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

% barsNylon.r = 0.05*ones(1,tData.nBar); % Radius of bars
% barsNylon.rho = 1150*ones(1,tData.nBar); % Density of bars (Nylon 6)
% barsNylon.nu = 0.46*ones(1,tData.nBar); % Poisson's ratio of bars (HDPE)
% barsNylon.E = 1e9*ones(1,tData.nBar); % Young's modulus of bars (HDPE)
% Point Masses
Mp = ones(1,tData.nPm); % All point masses initialised with a mass of 1

g = [0;0;0]; % Gravity

tData = tensegGenMat(tData,barsSoft,strings,Mp,g);
tDataRigid = tensegGenMat(tData,barsMetal,strings,Mp,g);

% The structure is in equilibrium in the initial position. Hence, there
% is no need to call the function tensegEq.

%% Simulation Inputs
tData.F = 1; % If 1, external forces are present in structure, not if 0.
tDataRigid.F = 1;
tData.damper = 100*ones(1,tData.nStr); % All strings initialised with dampers whose damping coefficient is 1. 
tDataRigid.damper = 100*ones(1,tData.nStr); % All strings initialised with dampers whose damping coefficient is 1. 
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

%% Compliance Matrix and Total Mass
Cred = inv(Kred); ecRed = eig(Cred);
CredRigid = inv(KredRigid); ecRedRigid = eig(CredRigid);
ex = 1:length(ecRed);
% figure();clf;
% scatter(ex, ecRed); hold on;
% scatter(ex, ecRedRigid);
% legend('soft','rigid');

mTotalSoft = 0; mTotalRigid = 0;
for k=1:tData.nBar
   mTotalSoft = mTotalSoft + tData.bars.listM(k);
   mTotalRigid = mTotalRigid + tDataRigid.bars.listM(k);
end
mTotalSoft
mTotalRigid
%% Deflection
% s = svd(Kred);    
% [U,~,V] = svd(Kred);
% sRigid = svd(KredRigid);
% [Urigid,~,Vrigid] = svd(KredRigid);
% 
% % ind = find(s > 1e-12);
% OneByK = zeros(ns-nC,ns-nC);
% OneByKRigid = zeros(ns-nC,ns-nC);
% for i = 1:length(s)
%    OneByK = OneByK + 1/s(i)*V(:,i)*(U(:,i).');
%    OneByKRigid = OneByKRigid + ...
%                    1/sRigid(i)*Vrigid(:,i)*(Urigid(:,i).');
% end
Fext = 100*ones(ns-nC,1);
Fext(1) = 0; Fext(3) = 0; %Fext(4) = 0;
% Fext(3) = 0; Fext(6) = 0; Fext(9) = 10; Fext(12) = 10;
% Fext(3) = 0; Fext(6) = 0; Fext(8) = 10; Fext(11) = 10;
% Fext(3) = 0; Fext(6) = 0; Fext(7) = 10; Fext(10) = 10;
deltx = Cred*Fext
deltxRigid = CredRigid*Fext

%%
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

options = odeset('RelTol',1e-10,'AbsTol',1e-10);%,'Refine',1);

% [simTime,tInt] = tensegSimTime(options,tEnd);

tData.Correction = 3; % Compressible Bar 
tDataRigid.Correction = 3; % Compressible Bar 
x0Sim = x0;
% x0Sim(9) = x0Sim(9) + 0.1*rand(1);
% x0Sim(12) = 0.1*rand(1);
%% Checking for H2
% extF = zeros(ns,length(simTime));
% for t=1:length(simTime)
%     extF(:,t) = [zeros(8,1);100*rand(1);zeros(3,1)];
% end
% tData.extF = extF;
% tDataRigid.extF = extF;
% [tFlex,yFlex] = tensegSim(x0Sim,simTime,tData,options);
% [tFlexRigid,yFlexRigid] = tensegSim(x0Sim,simTime,tDataRigid,options);
%% Linearization

[sysSS, Bf] = linSysCompDescriptor(x0, tData);
[sysSSRigid, BfRigid] = linSysCompDescriptor(x0, tDataRigid);

tData.sysSS = sysSS; tData.Bf = Bf;
tDataRigid.sysSS  = sysSSRigid; tDataRigid.Bf = BfRigid;

%% Simulate lqr for both soft and metal
% close all;
% tSim = (0:0.01:100)';
% x0Sim = zeros(2*ns,1);
% x0Sim(9) = 2*rand(1);
% x0Sim(12) = 2*rand(1);
% lqrLinSysComp(tData, tSim, x0Sim);
% lqrLinSysComp(tDataRigid, tSim, x0Sim);
%% Minreal
% Getting the stable part of the system
% sysSSRigid's subset of eigvec matrix not full rank
% buildStableNorm(sysSS)
% buildStableNorm(sysSSRigid)

minSys = minreal(sysSS);%,10);
minSysRigid = minreal(sysSSRigid);
% figure();
% pzmap(minSys); hold on;
% pzmap(minSysRigid);

%% Balred
% opts = balredOptions;
% pole(sysSS)
% opts = balredOptions('StateElimMethod','Truncate', 'Offset',1e-6,...
%         'AbsTol',1e-6, 'RelTol', 1e-6);
% bSys = balred(sysSS,18,opts);
% pole(bSys)

%% transforming sysSS using constraints
W2 = zeros(2*ns,2*(ns-nC));
W2(7,1) = 1; W2(9,2) = 1; W2(10,3) = 1; W2(12,4) = 1;
W2(19,5) = 1; W2(21,6) = 1; W2(22,7) = 1; W2(24,8) = 1;
Ared = W2.'*sysSS.A*W2; 
Bwred = W2.'*sysSS.B; Cred = W2.'*sysSS.C*W2; Dred = zeros(size(Bwred));
sysSSred = ss(Ared,Bwred,Cred,Dred);
AredRigid = W2.'*sysSSRigid.A*W2; 
BwredRigid = W2.'*sysSSRigid.B; CredRigid = W2.'*sysSSRigid.C*W2; 
DredRigid = zeros(size(Bwred));
sysSSredRigid = ss(AredRigid,BwredRigid,CredRigid,DredRigid);
eig(Ared)
eig(AredRigid)
%% stabsep
opt = stabsepOptions('AbsTol',1e-6,'Offset',1e-6);
[GS,GNS] = stabsep(sysSS,opt);
[GSRigid,GNS] = stabsep(sysSSRigid,opt);
pole(GS)
pole(GSRigid)
h2soft = norm(GS)
h2rigid = norm(GSRigid)
hinfsoft = norm(GS,Inf)
hinfRigid = norm(GSRigid,Inf)

Wc = gram(GS,'c'); Wo = gram(GS,'o');
WcRigid = gram(GSRigid,'c'); WoRigid = gram(GSRigid,'o');
h2gsoft = trace(GS.C*Wc*GS.C')
h2gsoft = trace(GS.B'*Wo*GS.B)
h2gRigid = trace(GSRigid.C*WcRigid*GSRigid.C')
h2gRigid =  trace(GSRigid.B'*WoRigid*GSRigid.B)
% figure();
% pzmap(GS); hold on;
% pzmap(GSRigid)
%%
% Amin = minSys.A;
% szAmin = size(Amin)
% rAmin = rank(Amin)
% 
% % [balSys,g] = balreal(sysSS)
% % rSys = balred(sysSS,szAmin(1))
% bodeplot(sysSS,minSys,'r--')
% Co = ctrb(minSys);
% unco = length(minSys.A) - rank(Co)


    

%% Checking Linear Model against Nonlinear
% x0delt = zeros(size(x0));%(1:end-1);
x0delt = zeros(2*ns,1);
delt(9) = 0.1*rand(1);
% x0delt = x0(1:end-1);% x0delt(9) = x0delt(9) + 0.1*rand(1);
T = 0:0.01:1;
for t = 1:length(T)
%     u(t,:) = [zeros(8,1);100*(sin(10*t)+sin(20*t)+cos(15*t)+cos(50*t))/4;zeros(3,1)];
    u(t,:) = [zeros(8,1);100*rand(1);zeros(3,1)];
end
% [tSim,ySim] = ode15s(@linearlagTensegrityDynamics_flex,T,x0delt,options,tData);
% [ysim,tsim,xsim] = lsim(sysSS,zeros(length(T),6),T,x0delt);
% [ysim,tsim,xsim] = lsim(sysSS,u,T,x0delt);
% [ysimRigid,tsim,xsim] = lsim(sysSSRigid,u,T,x0delt);
[ysim,tsim] = impulse(sysSS,T(end));
[ysimRigid,tsim] = impulse(sysSSRigid,T(end));
% % % plot(tSim,ySim(:,12));
% [tFlex,yFlex] = tensegSim(x0delt+x0,T,tData,options);
% for i=1:length(tFlex)
% %     ySim(i,:) = ySim(i,:) + x0(1:end)';
%     yFlex(i,:) = yFlex(i,:) - x0(1:end)'; % Comparing against eq
%     yFlexRigid(i,:) = yFlexRigid(i,:) - x0(1:end)'; % Comparing against eq
% end
% for t=1:length(tsim)
%     ysim(t,:) = ysim(t,:) - x0delt';
%     ysimRigid(t,:) = ysimRigid(t,:) - x0delt';
% end
%% Plotting 
% Plot Output Trajectories
plotMotion(tFlex, yFlex, tData);
plotMotion(tFlexRigid, yFlexRigid, tData);

% Plots to disturbance
% plotMotion(tsim, ysim, tData);
% plotMotion(tsim, ysimRigid, tDataRigid);

% Impulse response plots
% plotMotion(tsim, ysim(:,:,9), tData);
% plotMotion(tsim, ysimRigid(:,:,9), tDataRigid);



% Overlay Motion Plots
plotCompTbar_flex(tFlex,yFlex,tFlexRigid,yFlexRigid,tData);

% print(figure(3),'compTBarFlex_MotionNode3','-depsc');
% print(figure(4),'compTBarFlex_MotionNode4','-depsc');

% load('simSinRandomForce.mat')
% load('doubletNL.mat');