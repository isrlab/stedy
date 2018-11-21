function plot_configuration(q,tData,AZ,EL,axLims)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
% 
% Function to plot initial configuration of the tensegrity structure
% 
% INPUTS: [q, tData, AZ, EL]
% 
% q: vec(N) := initial nodal coordinates vectorized. 
% 
% AZ: Azimuth angle in degrees
%
% EL: Elevation angle in degrees
% 
% axLims: axis limits in figure window
% 

nNodes = length(q)/3;
figOptions = struct('pos',[450 40 500 500]);
figure(figOptions);hold on;

    % Draw the bars
    if(tData.nBar>0)
        for k=1:tData.nBar
            bk = tData.listX{k}*q;
            bkbar = tData.listXbar{k}*q;
            X = [(bkbar-bk/2) bkbar+bk/2];
            plot3(X(1,:), X(2,:), X(3,:),'k','Linewidth',2);
        end
    end
    
    % Draw the strings
    if(tData.nStr>0)
        for k=1:tData.nStr
            sk = tData.listY{k}*q;
            skbar = 0.5*abs(tData.listY{k})*q;
            X = [(skbar-sk/2) skbar+sk/2];
            plot3(X(1,:), X(2,:), X(3,:),'r','Linewidth',1);
        end
    end
    
    % Draw the point masses
    if(tData.nPm>0)
        for k=1:tData.nPm
            pk = tData.listP{k}*q;
            plot3(pk(1,:), pk(2,:), pk(3,:), 'k', 'Marker','square','MarkerSize',20);
        end
    end
    
    % Draw the nodes
    NN = reshape(q,3,nNodes);
    for k=1:nNodes
        if tData.fixedNodes(1,k)~=1
            plot3(NN(1,k),NN(2,k),NN(3,k),'ko','MarkerFaceColor','w')
        else
            plot3(NN(1,k),NN(2,k),NN(3,k),'ko','MarkerFaceColor','k')
        end
    end
    quiver3(zeros(3,1),zeros(3,1),zeros(3,1),[1;0;0],[0;1;0],[0;0;1],'Color','g')
    text(0.9,0,0,'x','fontsize',20,'color','b');
    text(0,0.9,0,'y','fontsize',20,'color','b');
    text(0,0,0.9,'z','fontsize',20,'color','b');
    
%     for i=1:tData.nNodes
%         txt{i} = num2str(i);
%         text(tData.N(1,i),tData.N(2,i),tData.N(3,i),txt{i},'fontsize',30);
%     end
%     axis equal;
%     axis off;
    view(AZ,EL);
end