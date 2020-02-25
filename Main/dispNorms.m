function dispNorms(tFlex, yFlex, tFlexRigid, yFlexRigid)
% Only for Tbar
%% Input signal
Ft = zeros(length(tFlex),12);
for i = 1:length(tFlex)
    t = tFlex(i);
    Ft(i,:) = 100*(t>=0 && t<2) - 100*(t>=2 && t<4) + 0*(t>=4);
end
enU = sqrt(trapz(tFlex,Ft(:,1).^2))
%% energy of output signal
% 3x and 3xdot
en3x = sqrt(trapz(tFlex,yFlex(:,7).^2));
en3xRigid = sqrt(trapz(tFlexRigid,yFlexRigid(:,7).^2));
en3xdot = sqrt(trapz(tFlex,yFlex(:,19).^2));
en3xdotRigid = sqrt(trapz(tFlexRigid,yFlexRigid(:,19).^2));

% 3z and 3zdot
en3z = sqrt(trapz(tFlex,yFlex(:,9).^2));
en3zRigid = sqrt(trapz(tFlexRigid,yFlexRigid(:,9).^2));
en3zdot = sqrt(trapz(tFlex,yFlex(:,21).^2));
en3zdotRigid = sqrt(trapz(tFlexRigid,yFlexRigid(:,21).^2));

% 4x and 4x dot
en4x = sqrt(trapz(tFlex,yFlex(:,10).^2));
en4xRigid = sqrt(trapz(tFlexRigid,yFlexRigid(:,10).^2));
en4xdot = sqrt(trapz(tFlex,yFlex(:,22).^2));
en4xdotRigid = sqrt(trapz(tFlexRigid,yFlexRigid(:,22).^2));

% 4x and 4x dot
en4z = sqrt(trapz(tFlex,yFlex(:,12).^2));
en4zRigid = sqrt(trapz(tFlexRigid,yFlexRigid(:,12).^2));
en4zdot = sqrt(trapz(tFlex,yFlex(:,24).^2));
en4zdotRigid = sqrt(trapz(tFlexRigid,yFlexRigid(:,24).^2));

enSoft = [en3x;en3z;en4x;en4z];
hinfsoft = max(enSoft)/enU
enRigid = [en3xRigid;en3zRigid;...
            en4xRigid;en4zRigid;];
hinfRigid = max(enRigid)/enU
%% max amplitude of output signal
max3x = max(yFlex(:,7));
max3xRigid = max(yFlexRigid(:,7));
max3xdot = max(yFlex(:,19));
max3xdotRigid = max(yFlexRigid(:,19));

max3z = max(yFlex(:,9));
max3zRigid = max(yFlexRigid(:,9));
max3zdot = max(yFlex(:,21));
max3zdotRigid = max(yFlexRigid(:,21));

max4x = max(yFlex(:,10));
max4xRigid = max(yFlexRigid(:,10));
max4xdot = max(yFlex(:,22));
max4xdotRigid = max(yFlexRigid(:,22));

max4z = max(yFlex(:,12));
max4zRigid = max(yFlexRigid(:,12));
max4zdot = max(yFlex(:,24));
max4zdotRigid = max(yFlexRigid(:,24));

maxSoft = [max3x;max3z;max4x;...
            max4z;];
h2soft = max(maxSoft)/enU        
maxRigid = [max3xRigid  ;max3zRigid ;...
            max4xRigid  ; max4zRigid ];
h2Rigid = max(maxRigid)/enU
