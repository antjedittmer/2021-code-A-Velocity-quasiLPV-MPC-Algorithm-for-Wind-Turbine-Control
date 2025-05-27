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
speedVec = [1,8,9,22];
figDirStr = 'figDir5';
[sysOut,gapCell] = compareLinearModels(speedVec,figDirStr); %
plotNormBodePlots(gapCell,speedVec,figDirStr);

%% Simulink simulations 
% Simulink models are compared with FAST (NREL) references.
% FAST simulation results obtained with FASTTool (Tu Delft) are provided as
% mat-files in dataIn folder. 

%Load data if available from previous simulation.
loadData = 0;
updateDDMdl1(0.75);

% Run Simulink models in closed loop w baseline controller( Torque controller
% k-omega-squared, Pitch controller: Gainscheduled Pi)
yAxCell = {'Wind (m/s)', 'GenTq (kNm)', 'Pitch (°)', 'RotSpd (rpm)',...
    'GenPwr (MW)','Twr_{FA} (m/s^2)', 'Twr_{SW} (m/s^2)'};

figNo = 2;
normStruct.Sweep = runCompareModels('Sweep',loadData,figNo,yAxCell,figDirStr);
figNo = figNo + 1;
normStruct.NTM18 = runCompareModels(18,loadData,figNo,yAxCell,figDirStr);
figNo = figNo + 1;
normStruct.EOG = runCompareModels('EOG',loadData,figNo,yAxCell,figDirStr);

save('normsGapCell','gapCell', 'normStruct');

%% Read out the quantative information
figNo = figNo + 1;
plotNormTimePlots(normStruct,figNo,figDirStr);

% Run models in closed loop with qLPV MPC
useFASTForComparison = 1;
figNo = figNo + 1;
runCompareCtrl('Sweep',loadData,figNo,useFASTForComparison,figDirStr);
figNo = figNo + 1;
runCompareCtrl('NTW18',loadData,figNo,useFASTForComparison,figDirStr);

