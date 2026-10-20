%% Script to generate pictures for paper
clc; clear; close all;

Simulink.data.dictionary.closeAll('-discard');

%% Set path to initialization script initWorkspace and run it
addpath(genpath('NonLinMdl'));
addpath(genpath('LinMdl'));

initWorkspace;

% Original models and qLPV-MPC of the paper (same results as the code of the
% paper): WECS model version 'paper' (getWecsModelVersion: no effective wind
% factor and Cq correction, FASTTool step reference; results in
% NonLinMdl/dataOut/paper) and the qLPV-MPC design 'paper' (qlmpcDesign:
% DTc = DT, no wind filter at the controller input, original formulation).
% The Cq scaling factors (updateDDMdl1) are set below.
% initWorkspace sets the default version 'FAST' again.
setModelVersionPaper('paper');

onlyDissPics = 1; % only generates diss plots, minus the filtered timeseries
allPlots = 0; % generate all plots, including norm plots}
filteredPlots = 1; % plots with filter for turbulent wind
% Rotor-effective wind filter (runCompareCtrl, pitch rate runs):
% 2: PI and qLMPC on the filtered plant Mdl_BianchiOL_WindFiltered, the
%    qLMPC scheduled with the filtered wind (setFilteredWindQLMPC)
% 1: as in the paper, PI on Mdl_BianchiOL_WindFiltered (turbulent wind
%    with the filter, stepped wind sweep with the filter bypassed, i.e. only
%    the projection factor kWindEffFilt), qLMPC on Mdl_BianchiOL
% 0: no filter, both controllers on the same plant Mdl_BianchiOL
% useFilteredCtrl: controller comparison (sweep and turbulent wind,
% cmpCtrlTimeDomain_All_RefSimulinkClosedLoop*_beta, results_tableSweep.tex,
% results_tableNTW18.tex, results_tableMerged.tex)
useFilteredCtrl = 0;
% addOldCtrl: additionally the controller comparison as before
% (useFiltered = 1, as in the paper), files with the suffix _old
addOldCtrl = 1;
% useFiltered: pitch rate runs below (PI on the filtered plant, qLMPC lines
% 4/8/13 deg/s without filter); additionally the PI without filter on
% Mdl_BianchiOL (cmpTimeDomain_All5PINewLim_NoFilt, statsTableMerged_NoFilt)
useFiltered = 1;
% Additional case (fourth case) with useFiltered = 2 for the turbulent wind:
% third column of the merged table (results_tableNTW18_windFilt.tex) and
% fourth line/column 'filt. (13)' of the qLMPC pitch rate figure
% and table (13 deg/s). See NonLinMdl/PlantControllerSetup.txt.
addFilteredCase = 1;

% Load data if available from previous simulation.
loadData = 1; % for runCompareModel (comparison Mdl1, Mdl2 and FAST)
loadDataCtrlTest = loadData; % for runCompareCtrl (PI and qLMPC based on Mdl2)
loadData1 = loadData; % for runCompareModels with qLMPC and PI at different rates


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


% Run Simulink models in closed loop w baseline controller( Torque controller
% k-omega-squared, Pitch controller: Gainscheduled Pi)
yAxCell = {'Wind (m/s)', 'GenTq (kNm)', 'Pitch (°)', 'RotSpd (rpm)',...
    'GenPwr (MW)','TwrAcc_{FA} (m/s^2)', 'TwrAcc_{SW} (m/s^2)'};

figNo = 2; % 4 m/s -> 
%updateDDMdl1(multRb,useTurn, x,scaleFactor)
updateDDMdl1(0.75, [],[], 0.96);
normStruct.Sweep = runCompareModels('Step',loadData,figNo,yAxCell,figDirStr,allPlots);
figNo = figNo + 2;
updateDDMdl1(0.75, [],[], 0.9475);
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
    runCompareModels(18,loadData,figNo,yAxCell,figDirStr,allPlots,ylimVal);

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

%% Comparison qLPV MPC vs. PI in closed loop
% 'SimulinkPaper': qLPV-MPC on Mdl_BianchiOL; PI on Mdl_BianchiOL
% (useFilteredCtrl = 0) or, as in the paper, on Mdl_BianchiOL_WindFiltered
% (_old files).
% Additional comparisons (runFASTComparison = 1) with the FAST-matched models
% (version 'FAST') and the qLPV-MPC design 'beta' of initWorkspace; both
% controllers run on the same plant model:
% - 'Simulink': simplified model Mdl_BianchiOL (uniform wind sweep)
% - 'SimFilt':  simplified model with rotor-effective wind filter
%               (Mdl_BianchiOL_WindFiltered) for the turbulent wind
% - 'FAST':     FAST v8.16 (NonLinMdl/FASTinputfiles; the 18 m/s case uses
%               the TurbSim file wind18_TurbKlasseA.bts)
% The qLPV-MPC settings are set in initWorkspace (qMPC, rMPC, pMPC, DTc, ...).
useFASTForComparison = 0; % compare with the PI on the same plant model
useTitle = [0,1]; % for thesis
ctrlCases = {'Sweep', 'SimulinkPaper'; 18, 'SimulinkPaper'};
% Pitch rate limit 13 deg/s for both controllers, set explicitly instead of
% using the values left in the data dictionaries by the rate loops below
setPitchRateLimitDD(13);
fileTableFilt = '';
if addFilteredCase % turbulent wind, both controllers on the filtered plant
    figNo = figNo + 1;
    runCompareCtrl(18, loadDataCtrlTest, figNo, useFASTForComparison, ...
        figDirStr, useTitle, 'SimulinkPaper', [], 2);
    fileTableFilt = 'results_tableNTW18_windFilt.tex';
end
% Comparison as before (as in the paper; the stepped wind sweep has no
% filtered PI run, so it is the PI with the filter bypassed): suffix _old
if addOldCtrl
    for idxC = 1:size(ctrlCases, 1)
        figNo = figNo + 1;
        runCompareCtrl(ctrlCases{idxC,1}, loadDataCtrlTest, figNo, useFASTForComparison, ...
            figDirStr, useTitle, ctrlCases{idxC,2}, [], 1, '_old');
    end
    merge_controller_tables('results_tableSweep_old.tex', 'results_tableNTW18_old.tex', ...
        'results_tableMerged_old.tex', fileTableFilt, 1);
end
for idxC = 1:size(ctrlCases, 1)
    figNo = figNo + 1;
    runCompareCtrl(ctrlCases{idxC,1}, loadDataCtrlTest, figNo, useFASTForComparison, ...
        figDirStr, useTitle, ctrlCases{idxC,2}, [], useFilteredCtrl);
end
merge_controller_tables('results_tableSweep.tex', 'results_tableNTW18.tex', ...
    'results_tableMerged.tex', fileTableFilt, useFilteredCtrl);

runFASTComparison = 0; % 1: additional comparisons on the FAST-matched models and FAST
if runFASTComparison
    setModelVersionPaper('FAST');
    ctrlCasesFAST = {'Sweep', 'Simulink'; 18, 'SimFilt'; 'Sweep', 'FAST'; 18, 'FAST'};
    for idxC = 1:size(ctrlCasesFAST, 1)
        figNo = figNo + 1;
        runCompareCtrl(ctrlCasesFAST{idxC,1}, loadDataCtrlTest, figNo, useFASTForComparison, ...
            figDirStr, useTitle, ctrlCasesFAST{idxC,2});
    end
    setModelVersionPaper('paper');
end

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
% Note: with loadData = 1 the qLPV-MPC results for the different rate
% limits (dataOut/OutTableTest_rate*_MPC.mat) are loaded if available; set
% loadData = 0 to recompute them with the current controller settings.
% The loops below leave the pitch rate limits in DD_Mdl1.sldd and
% DD_CtrlBaseline.sldd at the last value (13 deg/s).

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
controlObj = getEntry(dDataSectObj,'controlValue'); % read by MdlCtrl_qLPVMPC as controlValue.Pitch.Maxrate
controlValue = getValue(controlObj);

%--- Set vector which pitch rate constraints to be test
maxRateVector = [1,3,4,5,8,13];
outDataSimulationMat = 'OutDataWind18NTW.mat';
%dataDirOut = fullfile(fileparts(mfilename('fullpath')),'dataOut');
[~, ~, dataOutSub] = getWecsModelVersion;
dataDirOut = fullfile(fileparts(mfilename('fullpath')),'NonLinMdl','dataOut',dataOutSub); % same folder as the other simulation results

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

    % useFiltered = 2: qLPV-MPC on the filtered plant, scheduled with the
    % filtered wind (as runCompareCtrl; WECS model changed in memory only)
    strFilt = '';
    if useFiltered == 2
        oldValues = setFilteredWindQLMPC(simMdlname);
        if idx == 1, oldBaseVars = oldValues; end
        strFilt = '_windFilt';
    end

    if strcmp(outDataSimulationMat,'OutDataWind18NTW.mat')
        matFileOutTableTest1 = fullfile(dataDirOut,...
            sprintf('OutTableTest_rate%02d_MPC%s.mat',maxRate,strFilt));
    else
        matFileOutTableTest1 = fullfile(dataDirOut,...
            sprintf('OutTableTest_rate%02d_MPC_NTW16%s.mat',maxRate,strFilt));
    end

   [currentOutTableTestCell{idx}, tictoc_LPVMPCcell{idx},GenPwrRefCell{idx}] = ...
        getSimulationOutputTable(matFileOutTableTest1,loadData1,OutTable,simMdlname);
end
close_system(simMdlname, 0); % discard the WECS model change
if useFiltered == 2
    restoreBaseVars(oldBaseVars);
end

% Fourth case: 13 deg/s (last rate of the loop, still set in DD_Mdl1),
% qLPV-MPC on the filtered plant and scheduled with the filtered wind
OutTableMPCFilt13 = [];
if addFilteredCase && useFiltered ~= 2
    assert(maxRateVector(end) == 13, 'The fourth case needs 13 deg/s as the last rate.');
    load_system(simMdlname);
    oldBaseVars = setFilteredWindQLMPC(simMdlname);
    OutTableMPCFilt13 = getSimulationOutputTable(fullfile(dataDirOut, ...
        'OutTableTest_rate13_MPC_windFilt.mat'), loadData1, OutTable, simMdlname);
    close_system(simMdlname, 0);
    restoreBaseVars(oldBaseVars);
end


%% --- Get Data Directory value of pitch actuator rate
DDNameCtrl = 'DD_CtrlBaseline.sldd';
myDictionaryCtrlObj = Simulink.data.dictionary.open(DDNameCtrl);
dDataSectCtrlObj = getSection(myDictionaryCtrlObj,'Design Data');
controlCtrlObj = getEntry(dDataSectCtrlObj,'Control');
controlValue = getValue(controlCtrlObj);

currentOutTablePICell = cell(nRate,1); % cell with

simMdlname = 'test_SimulinkMdl2_Baseline'; % Bianchi model (includes blade state), same plant as the qLMPC runs

maxRateVector1 = 1:13;

% The PI test model uses the filtered plant Mdl_BianchiOL_WindFiltered; set
% the filter settings explicitly (as runCompareModels) instead of using the
% values left in the base workspace
if useFiltered >= 1
    filtValues = struct('step_thresshold', 1, 'kVeFilt', 0.25, 'kVeFilt_tau', 0.75);
    oldBaseVarsPI = struct;
    for fn = fieldnames(filtValues)'
        oldBaseVarsPI.(fn{1}) = evalin('base', fn{1});
    end
    restoreBaseVars(filtValues);
end

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
        getSimulationOutputTable(matFileOutTableTest1,loadData1,OutTable,simMdlname,...
        {'Wind', 'RotSpeed', 'GenPwr', 'GenTq', 'BlPitch1', 'NcIMUTAxs', 'NcIMUTAys','zetadotdot'}); % names are truncated to the number of logged signals
end
if useFiltered >= 1
    restoreBaseVars(oldBaseVarsPI);
end

%--- PI without the wind filter on the plant of the qLPV-MPC, Mdl_BianchiOL
% (no projection factor; as the controller comparison), at 4, 8 and
% 13 deg/s for cmpTimeDomain_All5PINewLim_NoFilt. The run at 13 deg/s is the
% same as the PI of the controller comparison (OutTableTest2NTW18_noFilt.mat).
maxRateVectorNoFilt = [4, 8, 13];
currentOutTablePINoFiltCell = cell(numel(maxRateVectorNoFilt), 1);
load_system(simMdlname);
blkPI = find_system(simMdlname, 'BlockType', 'ModelReference', 'Name', 'WECS Model');
set_param(blkPI{1}, 'ModelName', 'Mdl_BianchiOL'); % in memory only
for idx = 1:numel(maxRateVectorNoFilt)
    maxRate = maxRateVectorNoFilt(idx);
    controlValue = getValue(controlCtrlObj);
    controlValue.Pitch.Maxrate = maxRate;
    controlValue.Pitch.Minrate = - maxRate;
    setValue(controlCtrlObj,controlValue);
    saveChanges(myDictionaryCtrlObj)
    clear controlValue;

    matFileOutTableTest1 = fullfile(dataDirOut,...
        sprintf('OutTableTest_rate%02d_PI_noFilt.mat',maxRate));
    currentOutTablePINoFiltCell{idx} = ...
        getSimulationOutputTable(matFileOutTableTest1,loadData1,OutTable,simMdlname,...
        {'Wind', 'RotSpeed', 'GenPwr', 'GenTq', 'BlPitch1', 'NcIMUTAxs', 'NcIMUTAys','zetadotdot'});
end
close_system(simMdlname, 0); % discard the WECS model change

% qLPV-MPC on the filtered plant at 8 deg/s (result of the rate loop with
% useFiltered = 2), filtered case of statsTableMerged_NoFilt
fileMPCFilt8 = fullfile(dataDirOut, 'OutTableTest_rate08_MPC_windFilt.mat');
assert(isfile(fileMPCFilt8), 'Missing %s: run the qLPV-MPC rate loop once with useFiltered = 2.', fileMPCFilt8);
OutTableMPCFilt8 = getSimulationOutputTable(fileMPCFilt8, 1, OutTable, 'test_SimulinkMdl2_qLPVMPCbeta');

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

% input plotOutTable:OutTable,OutTableTest1,OutTableTest2,testCaseStr, testCaseCell, figNo2,figDir,strFig,timeEval)
timeEvalRate = [60, OutTable.Time(end)]; % evaluation window of the 18 m/s wind (60 to 660 s)
figDir = 'figDir';
testCaseStrMPC =  'qLMPC Rate Constraints (°/s)';
testCaseCell = {'standard (8)  ','high (13) ', 'low (4) ', 'filt. (13)'};
axPlotAllMPC = plotOutTable(tableForPlotMPC{1},tableForPlotMPC{2},tableForPlotMPC{3},...
   testCaseStrMPC,testCaseCell,[],figDir,[],timeEvalRate,OutTableMPCFilt13);
testCaseCell = testCaseCell(1:3); % PI: three cases

DDNameCtrl = 'DD_CtrlBaseline.sldd';
myDictionaryCtrlObj = Simulink.data.dictionary.open(DDNameCtrl);
dDataSectCtrlObj = getSection(myDictionaryCtrlObj,'Design Data');
controlCtrlObj = getEntry(dDataSectCtrlObj,'Control');
controlValue = getValue(controlCtrlObj);



% Plot qLmpc result
OutTable4 = currentOutTablePICell{maxRateVector1 == 4}; % PI results are indexed by maxRateVector1
OutTable4.Time = OutTable.Time(1:height(OutTable4));
testCaseStr =  'PI Rate Constraints (°/s)';

figNo2 = 200;
strFig = 'testConstrPI';
axPlotAll = plotOutTable(OutTable4,currentOutTablePICell{maxRateVector1 == 8},currentOutTablePICell{end},...
   testCaseStr,testCaseCell,figNo2,figDir,strFig,timeEvalRate);

% PI without filter (Mdl_BianchiOL) at 4, 8 and 13 deg/s; fourth case
% (purple): PI on the filtered plant at 8 deg/s. Statistics table
% figDir/statsTabletestConstrPI_NoFilt.tex
testCaseCellFilt8 = {'standard (8)  ','high (13) ', 'low (4) ', 'filt. (8)'};
OutTable4NoFilt = currentOutTablePINoFiltCell{maxRateVectorNoFilt == 4};
OutTable4NoFilt.Time = OutTable.Time(1:height(OutTable4NoFilt));
figNo3 = 201;
axPlotAllNoFilt = plotOutTable(OutTable4NoFilt, ...
    currentOutTablePINoFiltCell{maxRateVectorNoFilt == 8}, ...
    currentOutTablePINoFiltCell{maxRateVectorNoFilt == 13}, ...
    testCaseStr, testCaseCellFilt8, figNo3, figDir, 'testConstrPI_NoFilt', timeEvalRate, ...
    currentOutTablePICell{maxRateVector1 == 8});

% Set the axis for all axes: identical y limits in the qLMPC figure, the PI
% figure and the PI figure without filter
axCell = {axPlotAllMPC, axPlotAll, axPlotAllNoFilt};
for idx = 1:length(axPlotAllMPC)
    if all(cellfun(@(a) idx <= length(a) && isgraphics(a(idx)), axCell))
        lims = cell2mat(cellfun(@(a) get(a(idx),'YLim'), axCell, 'UniformOutput', false)');
        newLim = [min(lims(:,1)), max(lims(:,2))];
        cellfun(@(a) set(a(idx),'YLim',newLim), axCell);
    end
end

figDirConstr = 'figDirConstr1';
if ~isdir(figDirConstr)
    mkdir(figDirConstr);
end
figDir = figDirConstr;

strFig = 'MPCNewLim';
figMPC = ancestor(axPlotAllMPC(isgraphics(axPlotAllMPC)), 'figure');
if iscell(figMPC), figMPC = figMPC{1}; end
saveas(figMPC, fullfile(figDirConstr,['cmpTimeDomain_All5',strFig]), 'png');
saveas(figMPC, fullfile(figDirConstr, ['cmpTimeDomain_All5', strFig]), 'epsc');

strFig = 'PINewLim';
print(figNo2, fullfile(figDirConstr,['cmpTimeDomain_All5',strFig]), '-dpng');
print(figNo2, fullfile(figDirConstr, ['cmpTimeDomain_All5', strFig]), '-depsc');

strFig = 'PINewLim_NoFilt';
print(figNo3, fullfile(figDirConstr,['cmpTimeDomain_All5',strFig]), '-dpng');
print(figNo3, fullfile(figDirConstr, ['cmpTimeDomain_All5', strFig]), '-depsc');

% qLPV-MPC statistics with the filtered case at 8 deg/s
% (figDir/statsTabletestConstr_Filt8.tex; the figure is not used) and the
% merged table qLMPC and PI without filter: figDir/statsTableMerged_NoFilt.tex
figNo4 = 202;
plotOutTable(tableForPlotMPC{1},tableForPlotMPC{2},tableForPlotMPC{3},...
    testCaseStrMPC,testCaseCellFilt8,figNo4,'figDir','testConstr_Filt8',timeEvalRate,OutTableMPCFilt8);
merge_stats_tables(fullfile('figDir','statsTabletestConstr_Filt8.tex'), ...
    fullfile('figDir','statsTabletestConstrPI_NoFilt.tex'), ...
    fullfile('figDir','statsTableMerged_NoFilt.tex'), 'piNoFilt');

function setPitchRateLimitDD(maxRate)
% setPitchRateLimitDD sets the pitch rate limit maxRate (deg/s) of the PI
% (DD_CtrlBaseline.sldd, Control.Pitch) and of the qLPV-MPC
% (DD_Mdl1.sldd, controlValue.Pitch)
ddCell = {'DD_CtrlBaseline.sldd', 'Control'; 'DD_Mdl1.sldd', 'controlValue'};
for idxD = 1:size(ddCell, 1)
    dd = Simulink.data.dictionary.open(ddCell{idxD,1});
    entry = getEntry(getSection(dd, 'Design Data'), ddCell{idxD,2});
    val = getValue(entry);
    val.Pitch.Maxrate = maxRate;
    val.Pitch.Minrate = -maxRate;
    setValue(entry, val);
    saveChanges(dd);
end
end

function setModelVersionPaper(modelVersion)
% setModelVersionPaper sets the WECS model version ('paper' or 'FAST',
% getWecsModelVersion) with the corresponding qLPV-MPC design ('paper' or
% 'beta', qlmpcDesign) in the base workspace and updates the data dictionaries
assignin('base', 'wecsModelVersion', modelVersion);
if strcmp(modelVersion, 'paper')
    ctrlDesign = 'paper';
else
    ctrlDesign = 'beta';
end
designMPC = qlmpcDesign(ctrlDesign);
fn = fieldnames(designMPC);
for idxF = 1:numel(fn)
    assignin('base', fn{idxF}, designMPC.(fn{idxF}));
end
assignin('base', 'ctrlDesign', ctrlDesign);
assignin('base', 'designMPC', designMPC);
updateDDMdl1(0.75);
end
