clc; clear;

addpath(genpath('NonLinMdl'));
addpath(genpath('LinMdl'));

initWorkspace;

options = optimset('Display','iter','PlotFcns',@optimplotfval);

x0 = [2.1336    0.8631 1];


fun = @getInfluenceBladeParams;
[x,fval,exitflag,output] = fminsearch(fun,x0,options);

save('x','x');

% x =  [0.608380, 7128.534769];

% x =   [4.1814    0.7651]; % Optimized with wind speeds [1,22]

%x =  [2.1336    0.8631] % Optimized with wind speeds [1,22]


%x =  [2.1636    0.8885  -29.6035];