%% Script to generate pictures for paper
clc; clear; close all;

%% Set path to initialization script initWorkspace and run it
addpath(genpath('NonLinMdl'));
addpath(genpath('LinMdl'));

initWorkspace;

%% Matlab analysis: Bode plots linearized FAST models (reference)
% Bode plots of two linearized non-linear wind turbine models with linearized 
% FAST models. 

% Create Bode plots for comparison
compareLinearModels;

%% Simulink simulations 
% Simulink models are compared with FAST (NREL) references.
% FAST simulation results obtained with FASTTool (Tu Delft) are provided as
% mat-files in dataIn folder. 

%Load data if available from previous simulation.
loadData = 0;
updateDDMdl1(0.75);

% Run Simulink models in closed loop w baseline controller( Torque controller
% k-omega-squared, Pitch controller: Gainscheduled Pi)
figNo = 1;
normStruct.Sweep = runCompareModels('Sweep',loadData,figNo);
figNo = 4;
normStruct.NTM18 = runCompareModels(18,loadData,figNo);

%% Read out the quantative information

fieldsGapAll = fieldnames(normStruct); %get the field names
sortFieldsGap = fieldsGapAll; % Mdl1 (5 states) should be listed before Mdl2 (9 states)

sortFieldsSignals =  [{'GenTq'}  {'BlPitch1'}  {'RotSpeed'}    {'GenPwr'}   {'TwrAccFA'}    {'TwrAccSW'}];

% sysOutputname =[{'RotSpd \omega_r (-)'};{'TwrAcc_{fa} (-)'}; {'TwrAcc_{sw} (-)'}];
% sysInputname =[{'GenTq T_g (-)'}; {'BlPitch \beta_0 (-)'};{'Wind V_{\infty} (-)'}];

% NcIMUTAxs
windTickLabel = sortFieldsGap;


%% Set informations for both plots

lc = lines; %line colors
lc1 = [0*[0.5, 0.5,0.5];lc];
len = 6;
vec = 1:len;
idxVec = 1:len; % this is the number of input/output combination


%% Set values for first, norm plot
Mdl2Str = ' Rot,Twr,Gen+Bld';
legCell = ['                        Norm: \color{black}Mdl1: Rot+Twr ',...
    '\color[rgb]{',num2str(lc(1,:)),'}Mdl2:',Mdl2Str];

figure(100);
for idx = 1:6
    subplot(3,2,idx)
    aFieldname = sortFieldsSignals{idx};
    aFieldn = aFieldname;
    if strcmp(aFieldn,'TwrAccFA')
        aFieldn = 'NcIMUTAxs';
    elseif strcmp(aFieldn,'TwrAccSW')
         aFieldn = 'NcIMUTAys';
    end

    b = bar([normStruct.Sweep.(aFieldn); normStruct.NTM18.(aFieldn)]);

    for k = 1:length(b), b(k).FaceColor = lc1(k,:); end
    set(gca,'XTickLabel',sortFieldsGap);
    %set(gca,'YLim',[0,aMax]);

    if idx == 1
        title( legCell)
    end


    ylabel([aFieldname,' (-)']);

end

figFolder = 'figDir';
figStr = 'NormTime';
figFolderStr = fullfile(figFolder,figStr);
figFolderStrEps = figFolderStr;
print(figFolderStr, '-dpng');
print(figFolderStrEps, '-depsc');

% return;

% Run models in closed loop with qLPV MPC
figNo1 = 8;
runCompareCtrl('Sweep',loadData,figNo1);
figNo1 = 9;
runCompareCtrl('NTW18',loadData,figNo1);

