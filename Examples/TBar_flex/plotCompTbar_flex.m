function plotCompTbar_flex(t1,y1,t2,y2,tData)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
% 
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

for i=3:4
    figure(); clf;
    x_c = 3*(i-1) + 1;
    subplot(3,1,1);
    plot(t1,y1(:,x_c),t2,y2(:,x_c)); hold on;
%     plot(t1,abs(y1(:,x_c)-y2(:,x_c))); hold on;
    xlabel('Time (s)'); ylabel('x(t)'); 
    set(gca,'box','off');
    hold on;
    title(sprintf('Node %d', i));
    legend('Soft','Rigid');
    
    subplot(3,1,2);
    y_c = 3*(i-1) + 2;
    plot(t1,y1(:,y_c),t2,y2(:,y_c)); hold on;
%     plot(t1,abs(y1(:,y_c)-y2(:,y_c))); hold on;
    xlabel('Time (s)'); ylabel('y(t)'); 
    set(gca,'box','off');
    hold on;
    legend('Soft','Rigid');
%     legend('Rigid w/o Corr.','Compressible');
    
    subplot(3,1,3);
    z_c = 3*(i-1) + 3;
    plot(t1,y1(:,z_c),t2,y2(:,z_c)); hold on;
%     plot(t1,abs(y1(:,z_c)-y2(:,z_c))); hold on;
    xlabel('Time (s)'); ylabel('z(t)'); 
    set(gca,'box','off');
    hold on;
    legend('Soft','Rigid');
%     legend('Rigid w/o Corr.','Compressible');
end

end