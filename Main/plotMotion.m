function plotMotion(t,y,tData)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% Function to plot trajectories of all nodes in the tensegrity structure
% after solution has been obtained. 
%
% Each node's trajectory is plotted in a new figure with 3 subplots for
% x, y, and z projections. 
% 
% INPUTS: [t,y,tData]
% 
% t: Evaluation points, obtained from tensegSim.m
% 
% y: Solution matrix, obtained from tensegSim.m

for i=1:tData.nNodes
    figure(); clf;
    x_c = 3*(i-1) + 1;
    subplot(3,1,1);
    plot(t,y(:,x_c)); hold on;
    xlabel('Time (s)'); ylabel('x(t)'); 
    set(gca,'box','off');
    hold on;
    title(sprintf('Node %d', i));
    
    subplot(3,1,2);
    y_c = 3*(i-1) + 2;
    plot(t,y(:,y_c)); hold on;
    xlabel('Time (s)'); ylabel('y(t)'); 
    set(gca,'box','off');
    hold on;
    
    subplot(3,1,3);
    z_c = 3*(i-1) + 3;
    plot(t,y(:,z_c)); hold on;
    xlabel('Time (s)'); ylabel('z(t)'); 
    set(gca,'box','off');
    hold on;
end

end