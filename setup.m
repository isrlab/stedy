%% SETUP file to be run only the first time
% Copies ode45 and ode helper files from MATLAB's root directory to the current 
% working directory
% Adds all necessary functions to MATLAB path
% Edits ODE45 to include constraint correction function call

addpath([pwd filesep 'Main']);

cd Main
check = exist('ODE_helperFiles');
if(check == 0)
    mkdir ODE_helperFiles
end

srcLoc = [matlabroot filesep 'toolbox' filesep 'matlab' filesep 'funfun' filesep 'private'];
odeLoc = [matlabroot filesep 'toolbox' filesep 'matlab' filesep 'funfun' filesep 'ode45.m'];

copyfile(strcat(srcLoc),'ODE_helperFiles');
copyfile(strcat(odeLoc),'ode45m.m');
addpath([pwd filesep 'ODE_helperFiles']);

% Edit ode45m to add constraint correction function call and change
% function name
fileattrib ode45m.m +w;

fileID = fopen('ode45m.m','r');
formatSpec = '%c';
str = fscanf(fileID,formatSpec);
endStr ='      NNreset_f7 = false;';
strLocODE = strfind(str,endStr);
if(isempty(strLocODE))
    msg = 'Editing ode solver failed. Please edit ODE45m manually.';
    error(msg);
else 
    strCall = ['      ynew = ConstraintCorrection(ynew,odeArgs{1},tnew);' newline newline];
    str = insertBefore(str,endStr,strCall);
    
    % Change name
    endStr = '(ode,tspan';
    str = insertBefore(str,endStr,'m');
    fileID = fopen('ode45m.m','w+');
    formatSpec = '%s\n';
    fprintf(fileID,formatSpec,str);
    fclose(fileID);
end

pwd