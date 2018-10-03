function plotConstr(t,y,tData)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% Function to plot constraint violations over time.
% 
% INPUTS: [t,y,tData]
% 
% t: Evaluation points, obtained from tensegSim.m
% 
% y: Solution matrix, obtained from tensegSim.m
% 
% Each bar's length violation will be plotted in a new figure. 
%
% The last figure shall contain the extent of violation pertaining to the
% constraint on conservation of total energy.
% 

nn = numel(tData.N);

%% Plotting Errors in Bar Lengths over Time
bar_length = zeros(tData.nBar,length(t));
for j=1:tData.nBar
    for i=1:length(t)
        bar_length(j,i) = y(i,1:nn)*tData.NLin(j).Mat*y(i,1:nn)';
    end
end

for i=1:tData.nBar
    figure(); clf;
    semilogy(t,max(eps,sqrt((bar_length(i,1) - bar_length(i,:)).^2)));
    set(gca,'box','off');
    xlabel('Time (s)'); 
    ylabel(sprintf('Bar %d Length Error',i));
end

%% Plotting Errors in Total Energy over Time

T = 0*t;
V = 0*t; 
Vs = 0*t;
E_diff = 0*t;
E = 0*t;
work = y(:,end);
for i=1:length(t)
    q = y(i,1:3*tData.nNodes)';
    if(tData.nStr > 0)
        sk_norm = zeros(1,tData.nStr);
        for j=1:tData.nStr
            sk_norm(1,j) = norm(tData.listY{j}*q);
        end
        K = tData.K.*(sk_norm>tData.Lk); %% Check if strings are slack
        for j=1:tData.nStr
            Vs(i,1) = Vs(i,1)+0.5*K(j)*(sk_norm(j)-tData.Lk(j))^2;
        end
    else Vs(i,1) = 0;
    end
    T(i,1) = 1/2*y(i,3*tData.nNodes+1:end-1)*tData.M*y(i,3*tData.nNodes+1:end-1)';
    V(i,1) = -tData.G'*y(i,1:3*tData.nNodes)'+Vs(i,1);
    E(i,1) = T(i,1)+V(i,1);
    E_diff(i,1) = sqrt((-E(i)+work(i)+E(1))^2);
end

figure();clf;
semilogy(t,E_diff); xlabel('Time (s)'); ylabel('Error in Energy Conserved'); 
set(gca,'box','off');

end