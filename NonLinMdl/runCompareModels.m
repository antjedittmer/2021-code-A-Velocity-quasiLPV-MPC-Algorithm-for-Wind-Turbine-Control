function normStruct = runCompareModels(strWindType,loadData,figNo1,yAxCell)
% runCompareModels compares two Simulink models with FASTTool simulation 
% data generated with the baseline controller. 
%
% All inputs are optional:
% - strWindType: Two testcases: step sweep 4 to 25 ms and normal dist. with 
%   18 m/s mean (Default: Sweep)
% - loadData: load simulation output data if available instead of running
%   simulation (Default: 1)
% - figNo1: Number of figure (Default: 1)
% - yAxCell: Axes labels 

%% Handle optional inputs
% The default inputs are provided here.

% Two testcases: step sweep 4 to 25 ms and normal dist. 18 m/s mean
if ~nargin || isempty(strWindType) 
   strWindType = 'Sweep'; % 4,11,18,
end

if nargin < 2 || isempty(loadData) 
    loadData = 0; %load simulation output data if available;
end

if nargin < 3 || isempty(figNo1) 
    figNo1 = 1;
end
figNo2 = figNo1 + 1;

if nargin < 4 || isempty(yAxCell) % Axes labels for figure
    yAxCell = {'Wind V (m/s)', 'GenTq T_g (kNm)', 'Pitch \beta (°)', 'RotSpd \omega_r (rpm)',...
    'GenPwr P_g (MW)','Twr_{FA} y_t (m/s^2)', 'Twr_{SW} x_t (m/s^2)'};
end

%% Initialize path and files names
% The path to figure and input directory is set. The name of the Simulink
% models for Mdl1 and Mdl2 are provided. The name of the mat file with the
% FAST reference data is given.

% Set path to inputdata directory
workDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(workDir);
dataInDir = fullfile( mainDir,'dataIn');
addpath(dataInDir);

% Set path to figure directory
figDir = fullfile(mainDir,'figDir');
if ~isfolder(figDir)
    mkdir(figDir)
end

% Set path to data output directory
dataDirOut = fullfile(workDir, 'dataOut');
if ~isfolder(dataDirOut)
    mkdir(dataDirOut)
end

% Provide names of Simulink models to be run
simMdlname1 = 'test_SimulinkMdl1_Baseline'; 
simMdlname2 = 'test_SimulinkMdl2_Baseline'; 

% Provide names of FAST simulation data to be loaded
if strcmp(strWindType,'Sweep') == 1 % sweep from 4 to 25 in steps  
    outDataSimulationMat = 'OutDataSweep.mat'; %'OutDataStep.mat'; %
    strFig = '';
    testCaseStr = 'Wind Sweep';
else  % wind with average 18 m/s
    outDataSimulationMat = sprintf('OutDataWind%02dNTW.mat',strWindType);
    strFig = sprintf('NTW%02d',strWindType); %'NTW18';
    testCaseStr = sprintf('Wind, mean %02d m/s',strWindType);
end

%% Load data from FAST run
% Load FASTtool simulation data from dataIn folder.

load(outDataSimulationMat ,'OutTable');
if OutTable.Time(end) >1000
    idxT = OutTable.Time >= 100;
    OutTable = OutTable(idxT,:);
    OutTable.Time =  OutTable.Time - OutTable.Time(1);
end

% Assign data from FAST to variables
idxTime = 1: height(OutTable);
OutTable = OutTable(idxTime,:);
OutTable.BldPitch1 = OutTable.BlPitch1;
OutTable.Torque = OutTable.GenTq;

% Calculate and assign wind amplitude from comonents
idxWind = contains(OutTable.Properties.VariableNames, 'Wind1VelX');
vectWind = OutTable{:,idxWind};
vectAmpWind = sqrt(sum((vectWind.^2),2));

%% Run Simulink simulations or load data

% Names of simulation output in order of Simulink bus
% wind [m/s],Rotor Speed [rad/s],Generator Power [kW],Generator Torque  [Nm],Pitch [rad],xdotdotfa [m/s2],xdotdotsw [m/s2]
varnames = {'Wind', 'RotSpeed', 'GenPwr', 'GenTq', 'BlPitch1', 'NcIMUTAxs', 'NcIMUTAys','zetadotdot'}; 

% Run simulation for two models or load mat files if available
matFileOutTableTest1 = fullfile(dataDirOut,['OutTableTest1',strFig,'.mat']);
OutTableTest1 = getSimulationOutputTable(matFileOutTableTest1,loadData,OutTable,simMdlname1,varnames);

matFileOutTableTest2 = fullfile(dataDirOut,['OutTableTest2',strFig,'.mat']);
OutTableTest2 = getSimulationOutputTable(matFileOutTableTest2,loadData,OutTable,simMdlname2,varnames);

%% Create output plots (for use in power point)
% Two plots are created: 1st plot shows wind, rotor speed, 

% Create output plot time reference
idxTime = 1: min([length(idxTime), height(OutTableTest1),height(OutTableTest2)]);
time = OutTable.Time(idxTime);
idxTime = time > 60;
time = time(idxTime);

% Change line width
defaultLineWidth = get(groot,'defaultLineLineWidth');
set(groot,'defaultLineLineWidth',0.75);
cmpOnly1 = 2;


% Get color map for 'title legends'
cl = lines;
if cmpOnly1 <= 1
    titleStr = [testCaseStr,': {\color[rgb]{',num2str(cl(1,:)),'} Simulink Model, ',...
        '\color[rgb]{',num2str(cl(2,:)),'}FAST Model} '];
    titleStr1 = [testCaseStr,' Inputs: {\color[rgb]{',num2str(cl(1,:)),'} Simulink Model, ',...
        '\color[rgb]{',num2str(cl(2,:)),'}FAST Model} '];  
     titleStr2 = [testCaseStr,' Signals: {\color[rgb]{',num2str(cl(1,:)),'} Simulink Model, ',...
        '\color[rgb]{',num2str(cl(2,:)),'}FAST Model} '];
    
else
    titleStr = [testCaseStr,': Mdl1: Rot+Twr {\color[rgb]{',num2str(cl(1,:)),'}Mdl2: Rot,Gen,Twr+Bld ',...
        '\color[rgb]{',num2str(cl(2,:)),'}FAST} '];
    titleStr1 = titleStr;
    titleStr2 = titleStr;
end

% 1st figure
figure(figNo1);
axPlot(1) = subplot(3,1,1);
plot(time,vectAmpWind(idxTime),time,OutTableTest1.Wind(idxTime),'--'); 
axis tight; grid on;
ylabel('wind [m/s]')
title(titleStr1);

axPlot(2) = subplot(3,1,2);
if cmpOnly1 <= 1
    plot(time,OutTableTest2.GenTq(idxTime)/10^3,time,OutTable.GenTq(idxTime),'--');   
else
    plot(time,OutTableTest2.GenTq(idxTime)/10^3,time,OutTable.GenTq(idxTime),time,OutTableTest1.GenTq(idxTime)/10^3,'k--');
end
axis tight; grid on;
ylabel('GenTq T_g [kNm]') 

axPlot(3) = subplot(3,1,3);
if cmpOnly1 <= 1
    plot(time,OutTableTest2.BlPitch1(idxTime),time,OutTable.BlPitch1(idxTime),'--');
else
    plot(time,OutTableTest2.BlPitch1(idxTime),time,OutTable.BlPitch1(idxTime),time,OutTableTest1.BlPitch1(idxTime),'k--');
end
axis tight; grid on;
ylabel('BldPitch1 \beta [deg]')
xlabel('Time [s]')
linkaxes(axPlot,'x');
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)

if cmpOnly1 <= 1
    nameFig =  sprintf('cmpTimeDomain_WindIn%d%s',cmpOnly1,strFig);
else
    nameFig = ['cmpTimeDomain_WindIn',strFig];
end
set(gcf,'Name',['cmpTimeDomain_Wind',strFig])

print(fullfile(figDir,nameFig), '-dpng');

% 2nd figure
figure(figNo2)
axPlot2(1) = subplot(3,1,1);
if cmpOnly1 <= 1
    plot(time,OutTableTest2.RotSpeed(idxTime), time,OutTable.RotSpeed(idxTime),'--');
else
    plot(time,OutTableTest2.RotSpeed(idxTime), time,OutTable.RotSpeed(idxTime),time,OutTableTest1.RotSpeed(idxTime),'k--');
end
axis tight; grid on;
ylabel('RotSpd \omega_r [rpm]')
title(titleStr2);

axPlot2(2) = subplot(3,1,2);
if cmpOnly1 <= 1
    plot(time,OutTableTest2.NcIMUTAxs(idxTime),time,OutTable.NcIMUTAxs(idxTime),'--');
else
    plot(time,OutTableTest2.NcIMUTAxs(idxTime),time,OutTable.NcIMUTAxs(idxTime),time,OutTableTest1.NcIMUTAxs(idxTime),'k--');
end
axis tight; grid on;
ylabel('Twr_{FA} y_t[m/s^2]')

axPlot2(3) = subplot(3,1,3);
if cmpOnly1 <= 1
    plot(time,OutTableTest2.GenPwr(idxTime)/1000,time,OutTable.GenPwr(idxTime)/1000,'--');
else
    plot(time,OutTableTest2.GenPwr(idxTime)/1000,time,OutTable.GenPwr(idxTime)/1000,time,OutTableTest1.GenPwr(idxTime)/1000,'k--');
end
axis tight; grid on;
ylabel('GenPwr P_g [MW]')
xlabel('Time [s]')
set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
linkaxes(axPlot2,'x');

if cmpOnly1 <= 1
    nameFig = sprintf('cmpTimeDomain_WindSig1%d%s',cmpOnly1,strFig);
else
    nameFig = ['cmpTimeDomain_WindSig',strFig];
end

set(gcf,'Name',nameFig)

print(fullfile(figDir,nameFig), '-dpng');


%% L2 norm to quantify difference
varnamesnorm = varnames(2:7);
normStruct = struct;

for idx = 1: length(varnamesnorm)

    aVarname = varnamesnorm{idx};
    refVec = OutTable.(aVarname)(idxTime);
    if strcmp(aVarname,'GenTq')
        refVec = refVec* 1000;
    end

    normRef = norm(refVec);
    normStruct.(aVarname) = [norm(OutTableTest1.(aVarname)(idxTime) - refVec)/normRef,...
        norm(OutTableTest2.(aVarname)(idxTime) - refVec)/normRef];
end


%% Figure for paper and dissertation
figure(figNo2+1)
nAx = 7;
tiledlayout(nAx,1,'TileSpacing','Compact');

axPlotAll(1) = nexttile; %subplot(nAx,1,1);
plot(time,vectAmpWind(idxTime),time,OutTableTest1.Wind(idxTime),'k--'); 
axis tight; grid on;
ylabel(yAxCell{1}); %'wind [m/s]')
title(titleStr);
  
axPlotAll(2) = nexttile; %axPlotAll(2) = subplot(nAx,1,2);
plot(time,OutTableTest2.GenTq(idxTime)/10^3,time,OutTable.GenTq(idxTime),time,OutTableTest1.GenTq(idxTime)/10^3,'k--');
axis tight;
posAxis = axis;
axis([posAxis(1:2), min(43,posAxis(3)), 44]);
ylabel(yAxCell{3}); %'T_g [kNm]')
grid on; % axis tight;
%grid on;
ylabel(yAxCell{2}); %'GenTq T_g [Nm]') 

axPlotAll(3) = nexttile; %axPlotAll(3) = subplot(nAx,1,3);
plot(time,OutTableTest2.BlPitch1(idxTime),time,OutTable.BlPitch1(idxTime),time,OutTableTest1.BlPitch1(idxTime),'k--'); 
axis tight; grid on;
ylabel(yAxCell{3}); %'Pitch \beta [°]')

axPlotAll(4) = nexttile; %axPlotAll(4) = subplot(nAx,1,4);
plot(time,OutTableTest2.RotSpeed(idxTime), time,OutTable.RotSpeed(idxTime),time,OutTableTest1.RotSpeed(idxTime),'k--');
axis tight; grid on;
ylabel(yAxCell{4}); %'RotSpd \omega_r [rpm]')

axPlotAll(5) = nexttile; %axPlotAll(5) = subplot(nAx,1,5);
plot(time,OutTableTest2.NcIMUTAxs(idxTime),time,OutTable.NcIMUTAxs(idxTime),time,OutTableTest1.NcIMUTAxs(idxTime),'k--');
axis tight; grid on;
ylabel(yAxCell{6}); %'Twr_{FA} y_t[m/s^2]')

axPlotAll(6) = nexttile; %axPlotAll(6) = subplot(nAx,1,6);
plot(time,OutTableTest2.NcIMUTAys(idxTime),time,OutTable.NcIMUTAys(idxTime),time,OutTableTest1.NcIMUTAys(idxTime),'k--');
axis tight; grid on;
ylabel(yAxCell{7}); %'Twr_{SW} x_t[m/s^2]')
set(axPlotAll(6),'YLim', get(axPlotAll(5) ,'YLim'))

axPlotAll(7) = nexttile; %axPlotAll(7) = subplot(nAx,1,7);
plot(time,OutTableTest2.GenPwr(idxTime)/1000,time,OutTable.GenPwr(idxTime)/1000,time,OutTableTest1.GenPwr(idxTime)/1000,'k--');
axis tight; grid on;
ylabel(yAxCell{5}); %'GenPwr P_g [MW]')
xlabel('Time (s)')
linkaxes(axPlotAll,'x');

set(gcf,'Name',['cmpTimeDomain_All',strFig])
posDefault = get(0,'DefaultFigurePosition');
set(gcf, 'position', [posDefault(1),posDefault(2) - posDefault(4)*0.7,posDefault(3),posDefault(4)*2.1]);

set(findall(gcf,'-property','FontSize'),'FontSize',11.5)
set(findall(gcf,'-property','LineWidth'),'LineWidth',0.75)

%set(groot,'defaultLineLineWidth',defaultLineWidth);

print(fullfile(figDir,['cmpTimeDomain_All',strFig]), '-dpng');
print(fullfile(figDir,['cmpTimeDomain_All',strFig]), '-depsc');

