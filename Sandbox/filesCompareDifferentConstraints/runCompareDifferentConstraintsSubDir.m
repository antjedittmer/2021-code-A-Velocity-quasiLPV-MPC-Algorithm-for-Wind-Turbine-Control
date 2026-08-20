clc; clear; close all; bdclose all;
restoredefaultpath;

initWorkspaceSubDir;

matFileIn = 'OutDataWind18NTM.mat'; %'OutDataWind18NTW.mat'
%load(fullfile(dataInPath,matFileIn),'OutTable');
% dirMainProjects = fileparts(fileparts(fileparts(fileparts(fileparts(dataInPath)))));
% dirIPCtest = fullfile(dirMainProjects,'FASTTool','dataOut4');
% tmp = load(fullfile(dirIPCtest,'OutTable_IPC1_MaxRate4_Kp0dot00.mat'));
% OutTable = tmp.OutTable0;
% clear tmp;

%D:\ditt_aj\Projekte\FASTTool\dataOut4\OutTable_IPC1_MaxRate4_Kp0dot00.mat
%load(fullfile(dataInPath,'OutTable_TEF_off_IPC_off.mat'),'OutTable');

%addpath('subfunctionsMPC');

workDir = pwd;
mainDir = fileparts(workDir);

% Set path to figure directory
figDir = fullfile(mainDir,'figDirConstr');
if ~isfolder(figDir)
    mkdir(figDir)
end

%% Set path and parameter names
writeToExcel = 0;
% simMdlname = 'test_SimulinkMdl2_qLPVMPCbeta'; 
outDataSimulationMat = 'OutDataWind18NTM.mat'; %matFileIn;
loadData = 1;
onlyPlotCPCtimeseries = 0;

% Set path to data output directory
workDir = fileparts(mfilename('fullpath'));
dataDirOut = fullfile(workDir, 'dataOut');
if ~isfolder(dataDirOut)
    mkdir(dataDirOut)
end

%% Get Data Directory value of pitch actuator rate
DDName = 'DD_MdlCtrl_qLPVMPC.sldd';
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
controlObj = getEntry(dDataSectObj,'Control');
controlValue = getValue(controlObj);

%% Set vector which pitch rate constraints to be test
maxRateVector = 1:13;

%% MPC Loop over pitch rate constraint vector 
nRate = length(maxRateVector);
currentOutTableTestCell = cell(nRate,1); % cell with
tictoc_LPVMPCcell = cell(nRate,1);
GenPwrRefCell= cell(nRate,1);

simMdlname = 'test_SimulinkMdl2_qLPVMPCbeta'; %simMdlname = 'test_SimulinkMdl3LPV_MPC_a.slx';

for idx = 1:13
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

% Plot qLmpc result for different maximu rates
selR = [3,5,8,13];
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

figMPC = gcf;

%% Get Data Directory value of pitch actuator rate
DDNameCtrl = 'DD_CtrlBaseline.sldd';
myDictionaryCtrlObj = Simulink.data.dictionary.open(DDNameCtrl);
dDataSectCtrlObj = getSection(myDictionaryCtrlObj,'Design Data');
controlCtrlObj = getEntry(dDataSectCtrlObj,'Control');
controlValue = getValue(controlCtrlObj);

currentOutTablePICell = cell(nRate,1); % cell with

simMdlname = 'test_SimulinkMdl2_Baseline';
loadData1 = loadData;

for idx = 1: nRate
    maxRate = maxRateVector(idx);
    controlValue.Pitch.Maxrate = maxRate;
    controlValue.Pitch.Minrate = - maxRate;
    setValue(controlCtrlObj,controlValue);
    saveChanges(myDictionaryCtrlObj)
        
   
     if strcmp(outDataSimulationMat,'OutDataWind18NTW.mat')
        matFileOutTableTest1 = fullfile(dataDirOut,...
        sprintf('OutTableTest_rate%02d_PI.mat',maxRate));


     else
          matFileOutTableTest1 = fullfile(dataDirOut,...
        sprintf('OutTableTest_rate%02d_PI_NTW16.mat',maxRate));
    
    end
    
    currentOutTablePICell{idx} = ...
        getSimulationOutputTable(matFileOutTableTest1,loadData1,OutTable,simMdlname);
end

if writeToExcel == 1
    for idxE = 1:3
        spreadsheet = sprintf('MPC_CPC_%02d',selR(idxE+1));
        writetable(tableForPlotMPC{idxE},'NREL5MW_NTW18.xlsx','FileType','spreadsheet','Sheet',spreadsheet);
    end
    
    for idxE = 1:3
        spreadsheet = sprintf('PI_CPC_%02d',selR(idxE+1));
        writetable(currentOutTablePICell{maxRateVector == selR(idxE+1)},'NREL5MW_NTW18.xlsx','FileType','spreadsheet','Sheet',spreadsheet);
    end
end

% Plot qLmpc result
OutTable4 = currentOutTablePICell{maxRateVector == 4};
OutTable4.Time = OutTable.Time(1:height(OutTable4));
testCaseStr =  'PI Rate Constraints (°/s)';
figNo2 = 2;
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
saveas(figMPC, fullfile(figDirConstr,['cmpTimeDomain_All5',strFig]), 'png');
saveas(figMPC, fullfile(figDirConstr, ['cmpTimeDomain_All5', strFig]), 'epsc');

strFig = 'PINewLim';
print(gcf, fullfile(figDirConstr,['cmpTimeDomain_All5',strFig]), '-dpng');
print(gcf, fullfile(figDirConstr, ['cmpTimeDomain_All5', strFig]), '-depsc');

if onlyPlotCPCtimeseries == 1
    return
end

%% MPC limited to MPC2
tableForPlotMPC1{1} = currentOutTableTestCell{4};
tableForPlotMPC1{1}.Time = OutTable.Time(1:heightT);
tableForPlotMPC1{2} = currentOutTableTestCell{idx8};
tableForPlotMPC1{3} = currentOutTableTestCell{end};

testCaseStrMPC =  'qLMPC CPC Rate Constraints [deg/s]';
testCaseCell = {'standard (8) ','high (13) ', 'very low (2) '};
% axPlotAllMPC1 = plotOutTable(tableForPlotMPC1{1},tableForPlotMPC1{2},tableForPlotMPC1{3},...
%    testCaseStrMPC,testCaseCell,[],figDirConstr);
axPlotAllMPC1 = plotOutTable(OutTable,tableForPlotMPC1{2},tableForPlotMPC1{3},...
   testCaseStrMPC,testCaseCell,[],figDirConstr);


for idx = 1:length(axPlotAllMPC)
    if isgraphics(axPlotAll(idx))
        limPI = get(axPlotAll(idx),'YLim');
        limMPC = get(axPlotAllMPC1(idx),'YLim');
        newLim =[min(limPI(1),limMPC(1)), max(limPI(2),limMPC(2))];
        set(axPlotAllMPC1(idx),'YLim',newLim );
        set(axPlotAll(idx),'YLim',newLim );
    end
end
print(gcf, fullfile(figDir,['cmpTimeDomain_All5_2deg_s',strFig]), '-dpng');

%% IPC 
% load('OutTableIPC.mat','OutTableIPC');
% load('OutTableWeightsOrig.mat','OutTableCPC')

% tableForPlotIPC{1} = OutTableIPC;
% tableForPlotIPC{1}.Time = OutTable.Time(1:heightT);
% tableForPlotIPC{2} = currentOutTableTestCell{idx8};
% tableForPlotIPC{3} = currentOutTableTestCell{end};

% testCaseStrMPC =  'qLMPC CPC/IPC Rate Constraints [deg/s]';
% testCaseCell = {'CPC (8)','CPC (13)', 'IPC',};
% axPlotAllIPC1 = plotOutTable(tableForPlotIPC{1},tableForPlotIPC{2},tableForPlotIPC{3},...
%    testCaseStrMPC,testCaseCell);
% 
% for idx = 1:length(axPlotAllMPC)
%     if isgraphics(axPlotAll(idx))
%         limPI = get(axPlotAll(idx),'YLim');
%         limMPC = limPI; %get(axPlotAllMPC(idx),'YLim'); ToDo: Why is this not working
%         newLim =[min(limPI(1),limMPC(1)), max(limPI(2),limMPC(2))];
%         set(axPlotAllIPC1(idx),'YLim',newLim );
%         set(axPlotAll(idx),'YLim',newLim );
%     end
% end


%% Extract signals from table
idxT = 60/DT:height(currentOutTablePICell{1});
matrixPitchPI = cell2mat(cellfun(@(x) x.BlPitch1(idxT)', currentOutTablePICell, 'UniformOutput',false));
matrixPitchQLMPC = cell2mat(cellfun(@(x) x.BlPitch1(idxT)', currentOutTableTestCell, 'UniformOutput',false));
matrixTwrPI = cell2mat(cellfun(@(x) x.NcIMUTAxs(idxT)', currentOutTablePICell, 'UniformOutput',false));
matrixTwrQLMPC = cell2mat(cellfun(@(x) x.NcIMUTAxs(idxT)', currentOutTableTestCell, 'UniformOutput',false));

diffPitchPI = diff(matrixPitchPI');
diffPitchQLMPC = diff(matrixPitchQLMPC');

stdDiffPitchPI = mean(abs(diff(matrixPitchPI')))*180/pi;
stdDiffPitchQLMPC = mean(abs(diff(matrixPitchQLMPC')))*180/pi;
stdTowerPI = mean(abs(matrixTwrPI'));
stdTowerQLMPC = mean(abs(matrixTwrQLMPC'));

figure(3); 
vec4to13 = 1:3; % 3:13;
vec2to13 = 1:3; % 3:13;
subplot(2,2,1);
plot(maxRateVector(vec4to13),stdDiffPitchPI(vec4to13),'*-',...
    maxRateVector(vec2to13),stdDiffPitchQLMPC(vec2to13),'o-');
axis tight; grid on;
title('Mean values different rate limits')
ylabel('mean(\Delta Pitch) [deg/s]');

subplot(2,2,3);
plot(maxRateVector(vec4to13),stdTowerPI(vec4to13),'*-',...
    maxRateVector(vec2to13),stdTowerQLMPC(vec2to13),'o-');

axis tight; grid on;
ylabel('mean(TwrFA) [m/s^2]');
xlabel('Max. Pitch rate [m/s^2]');

subplot(1,2,2);
plot(stdDiffPitchPI(vec4to13),stdTowerPI(vec4to13),'*-',...
stdDiffPitchQLMPC(vec2to13),stdTowerQLMPC(vec2to13),'o-');
axis tight; grid on;
title({'Pareto front:','mean(TwrFa) vs. mean(\Delta Pitch)'})
legend('PI','qLMPC');
ylabel('mean(\Delta Pitch) [deg/s]');
xlabel('mean(TwrFA) [m/s^2]');

%% Barplots 

donotuse1Barplot = 0;
noTitle = 0;
aTestCell{1} = currentOutTablePICell{maxRateVector == 8};
aTestCell{1}.Time = OutTable.Time(1:heightT);
aTestCell{2} = currentOutTablePICell{maxRateVector == 13};
aTestCell{3} = currentOutTableTestCell{maxRateVector == 5};
aTestCell{3}.Time = OutTable.Time(1:heightT);
aTestCell{4-donotuse1Barplot} = currentOutTableTestCell{maxRateVector == 8};
aTestCell{5-donotuse1Barplot} = currentOutTableTestCell{maxRateVector == 13};
if donotuse1Barplot == 1
    summaryTickLabel = {'PI, limit 8°/s', 'PI, limit 13°/s',...
        'MPC, limit 4°/s','MPC, limit 13°/s'};
else
    summaryTickLabel = {'PI, limit 4°/s', 'PI, limit 13°/s','MPC, 4°/s',...
        'MPC, limit 8°/s','MPC, limit 13°/s'};
end
[K_matrix,strK] = calculateEvalCriteria(aTestCell,summaryTickLabel,noTitle);
figStr = 'PI4_13MPC4_13';
print(gcf,[fullfile(figDir,'optCriteria_PI_MPC'),'_',figStr], '-dpng');

printTitle = 0;
if strcmp(outDataSimulationMat,'OutDataWind16NTW.mat')
    save('matSubPlotBar_4to13_16.mat', 'K_matrix','strK','summaryTickLabel','printTitle');
else
    save('matSubPlotBar_4to13.mat', 'K_matrix','strK','summaryTickLabel','printTitle');
end

aNewTestCell{1} = aTestCell{1};
summaryTickLabelNew{1} = summaryTickLabel{1};
for idx = 1: length(currentOutTableTestCell)
aNewTestCell{idx+1} = currentOutTableTestCell{idx};
summaryTickLabelNew{idx+1} = sprintf('MPC, limit %d°/s',idx)
end
[K_matrix_all,strK1] = calculateEvalCriteria(aNewTestCell,summaryTickLabelNew,noTitle);

figure; plot(K_matrix_all(6:end,3),K_matrix_all(6:end,5))

figure; plot(K_matrix(3:end,3),K_matrix(3:end,5),'*-')


%[K_matrix,strk] = calculateEvalCriteria(aTestCell(3:5),summaryTickLabel(3:5),noTitle);

%load('matSubPlotBar.mat', 'K_matrix','strK','summaryTickLabel','printTitle');

% print(gcf,[fullfile(figDir,'cmpCtrlSimulink_PI_MPC'),'_',figStr,], '-depsc');

% for idx = 1:13, stdGP(idx) = std(currentOutTableTestCell{idx}.GenPwr); end
% for idx = 1:13, maxP(idx) = max(diff(currentOutTableTestCell{idx}.BlPitch1)/DT); end
% 
% for idx = 1:13, stdGP_PI(idx) = std(currentOutTablePICell{idx}.GenPwr); end
% for idx = 1:13, maxP_PI(idx) = max(diff(currentOutTablePICell{idx}.BlPitch1)/DT); end
% 
% 
% figure; plot(maxP,stdGP,'b-*'); hold on; plot(maxP_PI,stdGP_PI,'r-*');


%% Next bar plot: Only two bars omit for now
OutTable8 = currentOutTablePICell{end};
OutTable8.Time = OutTable.Time(1:height(OutTable8));

aTestCell1{1} = OutTable8;
aTestCell1{2} = currentOutTableTestCell{end};

summaryTickLabel = {'PI, limit 13°/s', 'MPC, limit 13°/s'};
calculateEvalCriteria(aTestCell1,summaryTickLabel,noTitle)
figStr = 'PI13MPC13';
print(gcf,[fullfile(figDir,'optCriteria_PI_MPC'),'_',figStr], '-dpng');

aTestCell3 = aTestCell1;
aTestCell3{2} = currentOutTableTestCell{end};
aTestCell3{3} = currentOutTableTestCell{2};

summaryTickLabel = {'PI, limit 13°/s', 'MPC, limit 13°/s', 'MPC, limit 2°/s'};
calculateEvalCriteria(aTestCell3,summaryTickLabel)

%[K_matrix_MPC, strK_MPC]
K_matrix_MPC = calculateEvalCriteria(aTestCell3,summaryTickLabel);


% load('OutTableIPC.mat','OutTableIPC');
load('OutTableWeightsOrig.mat','OutTableCPC')
cellOfTable{1} = OutTableCPC;
cellOfTable{2} = OutTableIPC;

% load('cellOfTable');
summaryTickLabel = {'PI, limit 8°/s', 'IPC, limit 8°/s'};

%%ellOfTable{3} = OutTableIPCSimAn;
[K_matrix_IPC, strK_IPC] = calculateEvalCriteria(cellOfTable,summaryTickLabel);

K_matrix_MPC_IPC = [K_matrix_IPC; K_matrix_MPC(2:end,:)];
strK.power = sprintf('%2.3f, ', K_matrix_MPC_IPC(:,1));
%strK.StdGenSpeed = sprintf('%2.3f, ',K_stdGenSpeed);
strK.StdPowSpeed = sprintf('%2.3f, ',K_matrix_MPC_IPC(:,2)); %K_stdGenPwr);
strK.ActPwr = sprintf('%2.3f, ',K_matrix_MPC_IPC(:,3)); %K_ActPwr);
strK.DamageBlades = sprintf('%2.3f, ',K_matrix_MPC_IPC(:,4)); %K_damageBlades);
strK.DamageTower = sprintf('%2.3f, ',K_matrix_MPC_IPC(:,5)); %K_damageTower);

summaryTickLabel = {'PI, 8°/s', 'PI-IPC, 8°/s','MPC, 13°/s', 'MPC, 2°/s'};

figDir = 'figDir';
plotCriteriaBar(K_matrix_MPC_IPC,strK,summaryTickLabel)
print(gcf,[fullfile(figDir,'optCriteria_PI_MPC_IPC_notitle'),'_',figStr], '-dpng');

K_matrix_MPC_IPC(3:4,:) = NaN;

plotCriteriaBar(K_matrix_MPC_IPC,strK,summaryTickLabel)
print(gcf,[fullfile(figDir,'optCriteria_PI_MPC_IPC1_notitle'),'_',figStr], '-dpng');

aTest = [currentOutTablePICell(8),tableForPlotMPC];
aCellforTest = [{'PI, limit 8°/s'}     {'MPC, 4° s'}    {'MPC, limit 8°/s'}    {'MPC, limit 13°/s'}];


figDir = '';

%% Plot old results
load('cellOfTable','cellOfTable');

cellOfTable4{1} = cellOfTable{1};
cellOfTable4{2} = cellOfTable{3};
strCell = {'PI-CPC','PI-IPC'}; %,'MPC, 13deg/s','MPC, 1deg/s'};
[K1_matrix_IPC, strK1_IPC] = calculateEvalCriteria(cellOfTable4,strCell);

% x = load('OutTableTest_rate13_IPC.mat');
%load('OutTableTest_rate08_PI.mat','OutTableTest_rate08_PI');
load('OutTableTest_rate01_MPC_LPV.mat','OutTableTest_rate01_MPC');
load('OutTableTest_rate13_MPC_LPV.mat','OutTableTest_rate13_MPC');
cellOfTable5 = cellOfTable4;
cellOfTable5{3} = OutTableTest_rate13_MPC;
cellOfTable5{4} = OutTableTest_rate01_MPC;
summaryTickLabel = {'PI, 8°/s', 'PI-IPC, 8°/s','MPC, limit 13°/s', 'MPC, limit 1°/s'};
[K_matrix_MPC, strK_MPC] = calculateEvalCriteria(cellOfTable5,summaryTickLabel);

% Plot qLmpc result
OutTable4 = OutTableCPC;
OutTable4.Wind = OutTable4.Wind1VelX;
OutTableIPC.Wind = OutTableIPC.Wind1VelX;
testCaseStr =  'Control Rate Constraints [deg/s]';
testCaseCell = {'CPC(8) ','IPC(13) ', 'low (1) ',};
plotOutTable(OutTable4,OutTableIPC,currentOutTableTestCell{end},...
   testCaseStr,testCaseCell);


tmp = load('cellOfTable');
x = load('OutTableTest_rate13_IPC.mat');
load('OutTableTest_rate08_PI.mat','OutTableTest_rate08_PI');
load('OutTableTest_rate01_MPC.mat','OutTableTest_rate01_MPC');
load('OutTableTest_rate13_MPC.mat','OutTableTest_rate13_MPC');





