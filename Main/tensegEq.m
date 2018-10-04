function tData = tensegEq(x0,x,tData)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% Function to find the equilibrium force densities in the strings at t=0. 
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
nLC = size(tData.Lin.A,1); % Number of linear constraints
nNLC = numel(tData.NLin); % Number of non-linear constraints
nConstr = nLC + nNLC;

gradLin = tData.Lin.A; % linear constraint,

gradNLin = [];

for i=1:nNLC
    gradNLin = [gradNLin;x'*tData.NLin(i).Mat];
end

gradR = [gradLin;gradNLin];

%% Setting up linear programming problem
n = numel(x);
Y = cell2mat(tData.Y);
L = reshape(Y'*x, n, tData.nStr);

Aeq = [L -gradR'];

if(tData.F)
    run('extF_eq');
    beq = tData.G + extF;
else beq = tData.G;
end

LB = [tData.minforce*ones(tData.nStr,1);-Inf*ones(nConstr,1)];

UB = Inf*ones(tData.nStr + nConstr,1);
options = optimoptions('linprog','Display','iter');
[SigLamb,fval,exitflag] = linprog([],[],[],Aeq,beq,LB,UB,[],options);

if(exitflag ~= 1)
    msg = 'Optimization failed. Please check input parameters again.';
    error(msg);
end

sigma = SigLamb(1:tData.nStr);

textSprFD = sprintf('%f \n',sigma); 
fprintf('Force Densities in Strings at Equilibrium: \n %s \n',textSprFD)

%% Computing Total Energy in Structure at Equilibrium
Vs = 0; % Potential Energy in Strings
if(tData.nStr>0)
    S = tData.N*tData.Cs'; % Strings
    for k=1:tData.nStr        
        s = S(:,k);
        L = norm(s); % Initial length between string nodes
        K(k) = tData.strings.E(k)*pi*(tData.strings.r(k))^2/L + sigma(k); % Stiffness of string
        Lk(k) = L*(1-sigma(k)/K(k));
        if Lk(k) < L
            Vs = Vs+0.5*K(k)*(L-Lk(k))^2;
        end
    end 
    tData.Lk = Lk;
    tData.K = K;
end 
tData.Vs = Vs;
tData.sigmaEq = sigma; % Force densities in springs at equilibrium

end