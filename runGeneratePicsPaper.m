%% Script to generate pictures for paper
clc; clear; close all;

Simulink.data.dictionary.closeAll('-discard');

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
loadData = 0;
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

%% For tests with filtered wind speed
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
        close_system(simMdlCell{idxM}, 0); % discard temporary edits (no .autosave)
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
useFASTForComparison = 0; % compare with PI on the same simplified Simulink model (Mdl_BianchiOL)
useTitle = [0,1]; % for thesis
loadDataCtrlTest = loadData; % this can be changed here
figNo = figNo + 1;
runCompareCtrl('Sweep',loadDataCtrlTest,figNo,useFASTForComparison,figDirStr,useTitle);
figNo = figNo + 1;
runCompareCtrl(18,loadDataCtrlTest,figNo,useFASTForComparison,figDirStr,useTitle);

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


%% Run baseline PI algorithm with different rate limits

%-- define sldd to change
%% Get Data Directory value of pitch actuator rate
DDName = 'DD_Mdl1.sldd'; %'DD_MdlCtrl_qLPVMPC.sldd';
loadData = loadDataCtrlTest;
try
    myDictionaryObj = Simulink.data.dictionary.open(DDName);
    if myDictionaryObj.HasUnsavedChanges
        discardChanges(myDictionaryObj);
    end
    myDictionaryObj.close();
catch
    % Dictionary wasn't open, nothing to do
end
myDictionaryObj = Simulink.data.dictionary.open(DDName);
dDataSectObj = getSection(myDictionaryObj,'Design Data');
controlObj = getEntry(dDataSectObj,'controlValue');
controlValue = getValue(controlObj);

%--- Set vector which pitch rate constraints to be test
maxRateVector = [1,3,4,5,8,13];
outDataSimulationMat = 'OutDataWind18NTW.mat';
%dataDirOut = fullfile(fileparts(mfilename('fullpath')),'dataOut');
dataDirOut = fullfile(pwd,'dataOut');

if strcmp(outDataSimulationMat, 'OutDataWind18NTW.mat')
    load(fullfile('dataIn',outDataSimulationMat),'OutTable');
end

%--- MPC Loop over pitch rate constraint vector 
nRate = length(maxRateVector);
currentOutTableTestCell = cell(nRate,1); % cell with
tictoc_LPVMPCcell = cell(nRate,1);
GenPwrRefCell= cell(nRate,1);

simMdlname = 'test_SimulinkMdl2_qLPVMPCbeta'; %simMdlname = 'test_SimulinkMdl3LPV_MPC_a.slx';

for idx = 1:length(maxRateVector)
    maxRate = maxRateVector(idx);
    controlValue.Pitch.Maxrate = maxRate;
    controlValue.Pitch.Minrate = - maxRate;
    setValue(controlObj,controlValue);
    saveChanges(myDictionaryObj);

   % Force model to pick up new dictionary values
    if bdIsLoaded(simMdlname)
        close_system(simMdlname, 0);  % close without saving
    end
    load_system(simMdlname);
    
    if strcmp(outDataSimulationMat,'OutDataWind18NTW.mat')
        matFileOutTableTest1 = fullfile(dataDirOut,...
            sprintf('OutTableTest_rate%02d_MPC.mat',maxRate));
    else
        matFileOutTableTest1 = fullfile(dataDirOut,...
            sprintf('OutTableTest_rate%02d_MPC_NTW16.mat',maxRate));
    end
       
   [currentOutTableTestCell{idx}, tictoc_LPVMPCcell{idx},GenPwrRefCell{idx}] = ...
        getSimulationOutputTable(matFileOutTableTest1,loadData,OutTable,simMdlname);
end


%% --- Get Data Directory value of pitch actuator rate
DDNameCtrl = 'DD_CtrlBaseline.sldd';
myDictionaryCtrlObj = Simulink.data.dictionary.open(DDNameCtrl);
dDataSectCtrlObj = getSection(myDictionaryCtrlObj,'Design Data');
controlCtrlObj = getEntry(dDataSectCtrlObj,'Control');
controlValue = getValue(controlCtrlObj);

currentOutTablePICell = cell(nRate,1); % cell with

simMdlname = 'test_SimulinkMdl1_Baseline';
loadData1 = 1; % loadData;

maxRateVector1 = 1:13;

for idx = maxRateVector1
    maxRate = maxRateVector1(idx);
    controlValue = getValue(controlCtrlObj);
    controlValue.Pitch.Maxrate = maxRate;
    controlValue.Pitch.Minrate = - maxRate;
    setValue(controlCtrlObj,controlValue);
    saveChanges(myDictionaryCtrlObj)
    clear controlValue;

    if strcmp(outDataSimulationMat,'OutDataWind18NTW.mat')
        matFileOutTableTest1 = fullfile(dataDirOut,...
            sprintf('OutTableTest_rate%02d_PI_1.mat',maxRate));
    else
        matFileOutTableTest1 = fullfile(dataDirOut,...
            sprintf('OutTableTest_rate%02d_PI_NTW16_1.mat',maxRate));
    end

    currentOutTablePICell{idx} = ...
        getSimulationOutputTable(matFileOutTableTest1,loadData1,OutTable,simMdlname);
end

% Plot qLmpc result for different maximu rates
selR = [3,4,8,13]; %[1,4,8,13]; % 
idx3 = maxRateVector == selR(1);
idx4 = maxRateVector == selR(2);
idx8 = maxRateVector == selR(3);
idx13 = maxRateVector == selR(4);

heightT = height(currentOutTableTestCell{idx4});
tableForPlotMPC{1} = currentOutTableTestCell{idx4}; 
tableForPlotMPC{1}.Time = OutTable.Time(1:heightT); % Time assigned to this table
tableForPlotMPC{2} = currentOutTableTestCell{idx8};
tableForPlotMPC{3} = currentOutTableTestCell{end};

% input plotOutTable:OutTable,OutTableTest1,OutTableTest2,testCaseStr, testCaseCell, figNo2,figDir,strFig)
figDir = 'figDir';
testCaseStrMPC =  'qLMPC Rate Constraints (°/s)';
testCaseCell = {'standard (8)  ','high (13) ', 'low (4) '};
axPlotAllMPC = plotOutTable(tableForPlotMPC{1},tableForPlotMPC{2},tableForPlotMPC{3},...
   testCaseStrMPC,testCaseCell,[],figDir);

DDNameCtrl = 'DD_CtrlBaseline.sldd';
myDictionaryCtrlObj = Simulink.data.dictionary.open(DDNameCtrl);
dDataSectCtrlObj = getSection(myDictionaryCtrlObj,'Design Data');
controlCtrlObj = getEntry(dDataSectCtrlObj,'Control');
controlValue = getValue(controlCtrlObj);



% Plot qLmpc result
OutTable4 = currentOutTablePICell{maxRateVector == 4};
OutTable4.Time = OutTable.Time(1:height(OutTable4));
testCaseStr =  'PI Rate Constraints (°/s)';

figNo2 = 200;
strFig = 'testConstrPI';
axPlotAll = plotOutTable(OutTable4,currentOutTablePICell{maxRateVector == 8},currentOutTablePICell{end},...
   testCaseStr,testCaseCell,figNo2,figDir,strFig);

% Set the axis for all axes
% xPlotAllMPC = axPlotAllMPC(isgraphics(axPlotAll, 'axes'));
% axPlotAll= axPlotAll(isgraphics(axPlotAll, 'axes'));
for idx = 1:length(axPlotAllMPC)
    if isgraphics(axPlotAll(idx))
        limPI = get(axPlotAll(idx),'YLim');
        limMPC = get(axPlotAllMPC(idx),'YLim');
        newLim =[min(limPI(1),limMPC(1)), max(limPI(2),limMPC(2))];
        set(axPlotAllMPC(idx),'YLim',newLim );
        set(axPlotAll(idx),'YLim',newLim );
    end
end

figDirConstr = 'figDirConstr1';
if ~isdir(figDirConstr)
    mkdir(figDirConstr);
end
figDir = figDirConstr;

strFig = 'MPCNewLim';
saveas(fullfile(figDirConstr,['cmpTimeDomain_All5',strFig]), 'png');
saveas(figDir, fullfile(figDirConstr, ['cmpTimeDomain_All5', strFig]), 'epsc');

strFig = 'PINewLim';
print(figNo2, fullfile(figDirConstr,['cmpTimeDomain_All5',strFig]), '-dpng');
print(figNo2, fullfile(figDirConstr, ['cmpTimeDomain_All5', strFig]), '-depsc');

