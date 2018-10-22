function animateTenseg(t,y,tData,time_loc,filename,format,frameRate,AZ,EL,axLims)
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

% Function to animate the simulated dynamics of the tensegrity structure
%
% INPUTS: [t, y, tData, time_loc, filename, format, frameRate, AZ, EL, AX]
% 
% t: Evaluation points, obtained from tensegSim.m
% 
% y: Solution matrix, obtained from tensegSim.m
% 
% time_loc: Location of text displaying current time to accompany the animation
% 
% filename: File name, specified as a character vector or string scalar.
% 
% format: Video format, like 'Uncompressed AVI' or 'MPEG-4'
% 
% frameRate: Rate of video playback in frames per second, specified as a positive number.
% 
% AZ: Azimuth angle in degrees
% 
% EL: Elevation angle in degrees
% 
% axLims: axis limits in animated figure
% 
if ismac
    vidObj = VideoWriter(filename,format);
elseif isunix
    sprintf('MPEG-4 is not a valid profile on Linux systems.');
    format = 'Motion JPEG AVI';
    filename = strcat(filename,'.avi');
    vidObj = VideoWriter(filename,format);
else
    vidObj = VideoWriter(filename,format);
end


vidObj.FrameRate = frameRate;
open(vidObj);

nNodes = size(tData.N,2);

h = figure('pos',[450 40 1200 1200]);
set(gcf, 'Color','white')
for i=1:length(t)
    figure(h); clf; hold on;
    q = y(i,1:3*nNodes)';
    
    % Draw the bars
    for k=1:tData.nBar
        bk = tData.listX{k}*q;
        bkbar = tData.listXbar{k}*q;
        X = [(bkbar-bk/2) bkbar+bk/2];
        plot3(X(1,:), X(2,:), X(3,:),'k','Linewidth',2);
    end
    
    % Draw the strings
    for k=1:tData.nStr
        sk = tData.listY{k}*q;
        skbar = 0.5*abs(tData.listY{k})*q;
        X = [(skbar-sk/2) skbar+sk/2];
        plot3(X(1,:), X(2,:), X(3,:),'r','Linewidth',1);
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
    
    % Display time
    text(time_loc(1),time_loc(2),time_loc(3),sprintf('Time: %.2f',t(i)),'fontsize',15)
    
    axis equal;
    axis(axLims);
    axis off
    view(AZ,EL)
    
    % drawnow;
    F = getframe(gca);
    writeVideo(vidObj,im2frame(F.cdata));

end
close(vidObj);