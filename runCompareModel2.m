function normStruct = runCompareModel2(strWindType,loadData,figNo1,yAxCell,figDirStr)
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
   strWindType = 'EOG'; % 4,11,18,
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

if nargin <5 || isempty(figDirStr)
    figDirStr = 'figDir';
end


%% Initialize path and files names
% The path to figure and input directory is set. The name of the Simulink
% models for Mdl1 and Mdl2 are provided. The name of the mat file with the
% FAST reference data is given.

% Set path to inputdata directory
workDir = fileparts(mfilename('fullpath'));
%mainDir = fileparts(workDir);
dataInDir = fullfile(workDir,'dataIn');
% addpath(dataInDir);


% Set path to data output directory
dataDirOut = fullfile(workDir, 'dataOut');
if ~isfolder(dataDirOut)
    mkdir(dataDirOut)
end

% Provide names of Simulink models to be run
% simMdlname1 = 'test_SimulinkMdl1_Baseline'; 
simMdlname2 = 'test_SimulinkMdl2_Baseline'; 

% Provide names of FAST simulation data to be loaded
if strcmp(strWindType,'Sweep') == 1 % sweep from 4 to 25 in steps  
    outDataSimulationMat = 'OutDataSweep.mat'; %'OutDataStep.mat'; %
    strFig = '';
    testCaseStr = 'Wind Sweep';
elseif isa(strWindType,'double')  % wind with average 18 m/s
    outDataSimulationMat = sprintf('OutDataWind%02dNTW.mat',strWindType);
    strFig = sprintf('NTW%02d',strWindType); %'NTW18';
    testCaseStr = sprintf('Wind, mean %02d m/s',strWindType);
else
     % EOG16mpers
      outDataSimulationMat = 'EOG16mpers.mat'; %'OutDataStep.mat'; %
    strFig = 'EOG';
    testCaseStr = 'EOG16mpers';
end

%% Load data from FAST run
% Load FASTtool simulation data from dataIn folder.

load(fullfile(dataInDir,outDataSimulationMat),'OutTable');
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
% matFileOutTableTest1 = fullfile(dataDirOut,['OutTableTest1',strFig,'.mat']);
% OutTableTest1 = getSimulationOutputTable(matFileOutTableTest1,loadData,OutTable,simMdlname1,varnames);

matFileOutTableTest2 = fullfile(dataDirOut,['OutTableTest2',strFig,'.mat']);
OutTableTest2 = getSimulationOutputTable(matFileOutTableTest2,loadData,OutTable,simMdlname2,varnames);

%% Create output plots (for use in power point)
% Two plots are created: 1st plot shows wind, rotor speed, 

% Create output plot time reference
idxTime = 1: min([length(idxTime),height(OutTableTest2)]);
time = OutTable.Time(idxTime);
idxTime = time > 60;
time = time(idxTime);

% Change line width
%defaultLineWidth = get(groot,'defaultLineLineWidth');
% set(groot,'defaultLineLineWidth',0.75);
% cmpOnly1 = 2;


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
    normStruct.(aVarname) = [norm(OutTableTest2.(aVarname)(idxTime) - refVec)/normRef];
end
