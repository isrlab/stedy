function tData = tensegConstruct(N,Cb,Cs,fixedNodes)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% Function to construct the essential parameters that uniquely define the
% tensegrity structure. 
%
% INPUTS: [N, Cb, Cs, fixedNodes]
%
% N: Initial nodal configuration, provided by the user. 
%
% Cb: Connectivity matrix corresponding to bars
%
% Cs: Connectivity matrix corresponding to strings
%
% fixedNodes: Matrix of the same size as N that details which nodes are
% fixed in which directions
%
% OUTPUT: tData
% tData is the MATLAB structure containing all fields requisite and
% pertinent to the dynamics simulation. 
%
% NEW FIELDS CREATED: nNodes, nBar, nStr, nPm, Lpm
%
% nNodes: Number of nodes in the tensegrity structure
%
% nBar: Number of bars in the structure, derived from Cb
%
% nStr: Number of strings in the structure, derived from Cs
%
% nPm: Number of point masses in the structure, placed exclusively at those
% nodes which do not connect to any bars
%
% Lpm: Location matrix of point masses, generated from Cb and N.

nNodes = size(N,2);

% Constructing Location matrix for point masses
Lp = zeros(nNodes,1); % Vector to identify which nodes are not connected to any bars
if(isempty(Cb)) % If there are no bars at all
    Lp = ones(nNodes,1); % All nodes with point masses
else
    for i=1:nNodes
        if(Cb(:,i) == zeros(size(Cb,1),1)) % If no bar connected to corresponding node
            Lp(i) = 1; % Placing a point mass at corresponding node
        end
    end
end
nPm = nnz(Lp); % Number of point masses
nzLp = find(Lp); % Indices in Lp with non-zero elements
Lpm = zeros(nPm,nNodes); % Location matrix of point masses

for i=1:nPm
    Lpm(i,nzLp(i)) = 1; % Matrix to identify locations of point masses:
                        % Size - number of point masses x number of Nodes
                        % In each row, 1 at node containing
                        % point mass and 0 eveywhere else
end

tData.N = N;
tData.nNodes = nNodes;
tData.Cb = Cb;
tData.Cs = Cs;
tData.nBar = size(Cb,1);
tData.nStr = size(Cs,1);
tData.nPm = nPm;
tData.Lpm = Lpm;
tData.fixedNodes = fixedNodes;

end