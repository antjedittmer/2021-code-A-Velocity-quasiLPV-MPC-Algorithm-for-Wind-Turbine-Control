%% Script to generate pictures for paper
clc; clear; close all;


%% Set path to initialization script initWorkspace and run it
addpath(genpath('NonLinMdl'));
addpath(genpath('LinMdl'));

initWorkspace;
onlyDissPics = 1; % only generates diss plots, minus the filtered timeseries
allPlots = 0; % generate all plots, including norm plots
filteredPlots = 1; % plots with filter for turbulent wind


%% Matlab analysis: Bode plots linearized FAST models (reference)
% Bode plots of two linearized non-linear wind turbine models with linearized
% FAST models.

% Create Bode plots for comparison
speedVec = [1,8,9,22] ;
figDirStr = 'figDir';

useActuatorStates = 0;
figNoAdd = 0;
createBodePlots = 1;
plotVisible = 'on';
noOut = 3;
xWeights = [1,1,1];

[sysOut,gapCell] = compareLinearModels(speedVec,figDirStr,...
    useActuatorStates,figNoAdd,createBodePlots,plotVisible,noOut,xWeights,allPlots); %

if onlyDissPics ==  0
    plotNormBodePlots(gapCell,speedVec,figDirStr);
end

%% Simulink simulations
% Simulink models are compared with FAST (NREL) references.
% FAST simulation results obtained with FASTTool (Tu Delft) are provided as
% mat-files in dataIn folder.

%Load data if available from previous simulation.
loadData = 1;
updateDDMdl1(0.75);

% Run Simulink models in closed loop w baseline controller( Torque controller
% k-omega-squared, Pitch controller: Gainscheduled Pi)
yAxCell = {'Wind (m/s)', 'GenTq (kNm)', 'Pitch (°)', 'RotSpd (rpm)',...
    'GenPwr (MW)','TwrAcc_{FA} (m/s^2)', 'TwrAcc_{SW} (m/s^2)'};

figNo = 2; % 4 m/s -> 
normStruct.Sweep = runCompareModels('Step',loadData,figNo,yAxCell,figDirStr,allPlots);
figNo = figNo + 2;
normStruct.Sweep = runCompareModels('Sweep',loadData,figNo,yAxCell,figDirStr,allPlots);
figNo = figNo + 1;
normStruct.NTM18 = runCompareModels(18,loadData,figNo,yAxCell,figDirStr,allPlots);

ax = findall(gcf, 'Type', 'Axes');
for idxA = numel(ax):-1:1
    yl(idxA,:) = ylim(ax(idxA));
    % fprintf('Axes %d: [%.4f, %.4f]\n', idxA, yl(idxA,1), yl(idxA,2)); %for debugging
end
ylimVal = flipud(yl);

% For tests with filtered wind speed
if filteredPlots == 1
    simMdlname1 = 'test_SimulinkMdl1_Baseline';
    simMdlname2 = 'test_SimulinkMdl2_Baseline';
    simMdlCell = {simMdlname1,simMdlname2};

    for idxM = 1: length(simMdlCell)
        load_system(simMdlCell{idxM})
        block = find_system(simMdlCell{idxM}, 'BlockType', 'ModelReference', 'Name', 'WECS Model');
        refModel = get_param(block{1}, 'ModelName');
        if ~contains(refModel,'WindFiltered')
            set_param(block{1}, 'ModelName', [refModel, '_WindFiltered']);
        end
    end
    figNo = figNo + 1;
    runCompareModels(18,0,figNo,yAxCell,figDirStr,allPlots,ylimVal);

    % Remove the filters again
    for idxM = 1: length(simMdlCell)
        load_system(simMdlCell{idxM})
        block = find_system(simMdlCell{idxM}, 'BlockType', 'ModelReference', 'Name', 'WECS Model');
        refModel = get_param(block{1}, 'ModelName');
        set_param(block{1}, 'ModelName', strrep(refModel, '_WindFiltered',''));

    end
end

if onlyDissPics == 0
    figNo = figNo + 1;
    normStruct.EOG = runCompareModels('EOG',loadData,figNo,yAxCell,figDirStr);
    save('normsGapCell','gapCell', 'normStruct');
end


%% Read out the quantative information
% figNo = figNo + 1;
% plotNormTimePlots(normStruct,figNo,figDirStr);

% Run models in closed loop with qLPV MPC
useFASTForComparison = 1;
useTitle = [0,1]; % for thesis
loadDataCtrlTest = loadData; % this can be changed here
figNo = figNo + 1;
runCompareCtrl('Sweep',loadData,figNo,useFASTForComparison,figDirStr,useTitle);
figNo = figNo + 1;
runCompareCtrl(18,loadData,figNo,useFASTForComparison,figDirStr,useTitle);

%% For debugging

% varnames = {'Wind', 'RotSpeed', 'GenPwr', 'GenTq', 'BlPitch1', ...
%     'NcIMUTAxs', 'NcIMUTAys'};
% OutDataTable13 = array2table(OutDataTest,'VariableNames',varnames);
% time0 = OutTable.Time(1:length(OutDataTable.GenPwr));
% idxT = time0 >= 30;
% time1 = time0(idxT);
% figure; subplot(2,1,1); plot(time1,OutDataTable13.GenPwr(idxT));
% subplot(2,1,2); plot(time1,OutDataTable13.NcIMUTAxs(idxT));
% OutDataTable04= array2table(OutDataTest,'VariableNames',varnames);
% OutDataTable08= array2table(OutDataTest,'VariableNames',varnames);

