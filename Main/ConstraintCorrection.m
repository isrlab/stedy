function ynew = ConstraintCorrection(x,tData,t)
% This Source Code Form is subject to the terms of the Mozilla Public
% License, v. 2.0. If a copy of the MPL was not distributed with this
% file, You can obtain one at http://mozilla.org/MPL/2.0/.

% Function to minimise constraint violations occurring due to numerical
% integration. Called after a successful step from ODE45.m
%
% INPUTS: [x,tData,t]
%
% x: Solution vector at time t
%
% tData.Correction: 0 or 1. If 1, inclusive of total energy correction. 
%

    work = x(end);
    ns = (numel(x)-1)/2;
    q = x(1:ns); % Position vector
    qd = x(ns+1:end-1); % Velocity vector
    
    gradLin = tData.Lin.A; % linear constraint,
    nLC = size(tData.Lin.A,1); % Number of linear constraints
    nNLC = numel(tData.NLin); % Number of non-linear constraints
    nConstr = nLC + nNLC;
    
    gradNLin = [];   
    for ii=1:nNLC
        gradNLin = [gradNLin;2*q'*tData.NLin(ii).Mat];
    end
    gradR = [gradLin;gradNLin];
    
    II = eye(ns);
    X = zeros(nConstr,ns);
    for i=1:ns
        gradNLindq = zeros(nLC,ns);
        for ii=1:nNLC
            gradNLindq = [gradNLindq;2*tData.NLin(ii).Mat(i,:)];
        end
        X = X+gradNLindq*qd*II(i,:);
    end
    
    % Check to see if bars are present in structure
    if(nNLC > 0)
        R = [gradLin;0.5*gradNLin]*q-[tData.Lin.b;tData.bars.L0.^2];
    else
        R = gradLin*q - tData.Lin.b;
    end

ynew = x;

if(tData.Correction == 1) % With Energy Correction
    
    Vs = 0; % Potential energy stored in strings
    if(tData.nStr>0)
        for k=1:tData.nStr
            s = tData.listY{k}*q;
            L = norm(s); % Initial length between string nodes
            if tData.Lk(k)<L
                Vs = Vs+0.5*tData.K(k)*(L-tData.Lk(k))^2;
            end
        end 
    end

    Energy = 0.5*x(ns+1:end-1)'*tData.M*x(ns+1:end-1)-tData.G'*x(1:ns)+Vs; % Total energy in structure
    
    while (norm(R) > 10^-10 || abs(tData.Energy + work - Energy) > 10^-10 )   
        if tData.nStr>0 
            for i=1:tData.nStr
                sk_norm(1,i) = norm(tData.listY{i}*q);
            end
            K = tData.K.*(sk_norm>tData.Lk); %% Check if strings are slack
            E_q = -tData.G'+q'*cell2mat(tData.Y)*kron(K',eye(ns))-q'*cell2mat(tData.Y)*kron(K'.*tData.Lk'./sk_norm',eye(ns));
            E_dq = qd'*tData.M;
            E0 = 0.5*qd'*tData.M*qd-tData.G'*q+0.5*q'*cell2mat(tData.Y)*kron(K',eye(ns))*q-sum(K.*tData.Lk.*sk_norm)+0.5*sum(K.*tData.Lk.^2);
        else

            E_q = -tData.G';
            E_dq = qd'*tData.M;
            E0 = 0.5*qd'*tData.M*qd-tData.G'*q;
        end

        A = [gradR zeros(size(gradR)) zeros(size(gradR,1),1);X gradR zeros(size(gradR,1),1);E_q E_dq -1];
        ynew = A'/(A*A')*[-R;-gradR*qd;tData.Energy-E0+(work)] + x;
          
        %% To calculate updated constraint violations
        x = ynew;
        work = x(end);
        ns = (numel(x)-1)/2;
        q = x(1:ns);
        qd = x(ns+1:end-1);
        gradLin = tData.Lin.A; % linear constraint,
        nLC = size(tData.Lin.A,1);
        nNLC = numel(tData.NLin);
        nConstr = nLC + nNLC;

        gradNLin = [];   
        for ii=1:nNLC
            gradNLin = [gradNLin;2*q'*tData.NLin(ii).Mat];
        end
        gradR = [gradLin;gradNLin];

        II = eye(ns);
        X = zeros(nConstr,ns);
        for i=1:ns
            gradNLindq = zeros(nLC,ns);
            for ii=1:nNLC
                gradNLindq = [gradNLindq;2*tData.NLin(ii).Mat(i,:)];
            end
            X = X+gradNLindq*qd*II(i,:);
        end

        % added to be able to run cases with no bars
        if(nNLC > 0)
            R = [gradLin;0.5*gradNLin]*q-[tData.Lin.b;tData.bars.L0.^2];
        else
            R = gradLin*q - tData.Lin.b;
        end

        Vs = 0; % Potential energy stored in strings
        if(tData.nStr>0)
            for k=1:tData.nStr
                s = tData.listY{k}*q;
                L = norm(s); % Initial length between string nodes
                if tData.Lk(k)<L
                    Vs = Vs+0.5*tData.K(k)*(L-tData.Lk(k))^2;
                end
            end 
        end
        
        Energy = 0.5*x(ns+1:end-1)'*tData.M*x(ns+1:end-1)-tData.G'*x(1:ns)+Vs; % Total energy in structure
    end
    
else % No Energy Correction
    
    while (norm(R) > 10^-10)   
        A = [gradR zeros(size(gradR));X gradR]; % Without Energy correction
        ynew = A'/(A*A')*[-R;-gradR*qd] + x(1:end-1); % Without Energy correction
        ynew = [ynew;0]; % Without Energy correction
                
        %% To calculate updated constraint violations
        x = ynew;
        work = x(end);
        ns = (numel(x)-1)/2;
        q = x(1:ns);
        qd = x(ns+1:end-1);
        gradLin = tData.Lin.A; % linear constraint,
        nLC = size(tData.Lin.A,1);
        nNLC = numel(tData.NLin);
        nConstr = nLC + nNLC;

        gradNLin = [];   
        for ii=1:nNLC
            gradNLin = [gradNLin;2*q'*tData.NLin(ii).Mat];
        end
        gradR = [gradLin;gradNLin];

        II = eye(ns);
        X = zeros(nConstr,ns);
        for i=1:ns
            gradNLindq = zeros(nLC,ns);
            for ii=1:nNLC
                gradNLindq = [gradNLindq;2*tData.NLin(ii).Mat(i,:)];
            end
            X = X+gradNLindq*qd*II(i,:);
        end

        % added to be able to run cases with no bars
        if(nNLC > 0)
            R = [gradLin;0.5*gradNLin]*q-[tData.Lin.b;tData.bars.L0.^2];
        else
            R = gradLin*q - tData.Lin.b;
        end
    end
end