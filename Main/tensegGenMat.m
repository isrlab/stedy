function tData = tensegGenMat(tData,bars,strings,Mp,g)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% Function to generate matrices under the Lagrangian framework described in
% the theoretical paper. 
%
% INPUTS: [tData, bars, strings, Mp, g]
%
% bars: MATLAB structure with fields describing geometric and material
% properties of the bars in the structure. 
% 
% strings: MATLAB structure with fields describing geometric and material
% properties of the strings in the structure. 
%
% Mp: Vector containing masses of the point masses in the structure. 
%
% g: Gravity vector
%
% 
% NEW FIELDS CREATED in tData: [M, G, listX, listXbar, listP, Nlin,...
%                               ...Lin, Yk, Lk, K, Y, Vs, bars.L0]
%
% M: Mass matrix from Total Kinetic Energy
% 
% G: from Gravitational Potential Energy
% 
% listX: Cell with each element containing matrix Xk to help coordinatise
% bar bk
% 
% listXbar: Cell with each element containing matrix Xkbar to help
% coordinatise center of mass of bar bk
% 
% listP: Cell with each element containing matrix Pk to help coordinatise
% point of mass pk
%
% Nlin: MATLAB Structure array with matrix fields corresponding to non-linear
% constraints, i.e., bar lengths
%
% Lin: MATLAB structure array with matrix fields corresponding to linear
% constraints, i.e., fixed nodes.
%
% Yk: Cell with each element containing matrix Yk to help coordinatise
% string sk
% 
% Lk: Vecor containing rest lengths of all strings in the structure. 
%
% K: Vector containing stiffness of all strings in the structure. 
%
% Y: MATLAB structure array with matrix fields corresponding to Yk'*Yk to
% help compute potential energy stored in strings
% 
% Vs: Scalar field containing potential energy stored in strings. 
% 
% bars.L0: Inital fixed length of the bars in the structure. 
% 

tData.g = g;

TH = eye(tData.nBar); % Theta matrix of bars
PSI = eye(tData.nStr); % Psi matrix of strings
PHI = eye(tData.nPm); % Phi matrix of point masses

B = tData.N*tData.Cb'; % Bars
S = tData.N*tData.Cs'; % Strings
P = tData.N*tData.Lpm'; % Point masses

M = zeros(3*tData.nNodes,3*tData.nNodes);
G = zeros(3*tData.nNodes,1);

% Linear Constraints Aq = b
fixedA = diag(tData.fixedNodes(:)); 
ii = find(sum(fixedA,2) ~= 0);
Lin.A = fixedA(ii,:);
bb = tData.N(:);
Lin.b = bb(ii);

% Bars
for k=1:tData.nBar
    thk = TH(:,k); 
    Xk = kron(thk'*tData.Cb,eye(3));
    Xkbar = 0.5*kron(thk'*abs(tData.Cb),eye(3));
    b = B(:,k);
    L = norm(b);
    mk = bars.rho(k)*(pi*bars.r(k)^2)*L;
    Ik = (mk/12)*(L^2); 
    M = M + mk*Xkbar'*Xkbar + (Ik/L^2)*Xk'*Xk;
    G = G + mk*Xkbar'*tData.g; 
    NLin(k).Mat = Xk'*Xk; % Nonlinear constraints pertaining to bar lengths
    bars.L0(k,1) = L; % Initial lengths of bars
    tData.listX{k} = Xk;
    tData.listXbar{k} = Xkbar;
end

% Point Masses
if(tData.nPm>0) 
    for k=1:tData.nPm
        phiK = PHI(:,k);
        Pk = kron((phiK'*tData.Lpm),eye(3));
        p = P(:,k); 
        mp = Mp(k);
        M = M + mp*Pk'*Pk;
        G = G + mp*Pk'*g;
        tData.listP{k} = Pk;
    end
end
tData.bars = bars;
tData.M = M; % Matrix from Kinetic energy expression
tData.G = G; % Matrix from gravitational potential energy expression
tData.NLin = NLin; % Non-linear constraints
tData.Lin = Lin; % Linear constraints

% Strings
Vs = 0; % Potential energy stored in strings
if(tData.nStr>0)
    for k=1:tData.nStr
        psiK = PSI(:,k);
        Yk = kron((psiK'*tData.Cs),eye(3));
        s = S(:,k);
        L = norm(s); % Initial length between string nodes
        Lk(k) = strings.rLP(k)*L; % Rest Length of String (initial)
        K(k) = strings.E(k)*pi*(strings.r(k))^2/Lk(k); % Stiffness of string
        Y{k} = Yk'*Yk;
        if Lk(k)<L
            Vs = Vs+0.5*K(k)*(L-Lk(k))^2;
        end
        tData.listY{k} = Yk;
    end 
    tData.Y = Y;
    tData.Lk = Lk;
    tData.K = K;
end

tData.Vs = Vs;

tData.bars = bars;
tData.strings = strings;

end