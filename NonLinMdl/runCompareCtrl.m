function  OutTableTest2  = runCompareCtrl(strWindType,loadData,figNo1,useFASTForComparison,figDirStr,useTitle,plantMdl,ctrlOpts,useFiltered,strOut)
% runComparCtrl compares simulation results for baseline and qLPV MPC
% controller in closed loop with the Simulink WECS model.
%
% All inputs are optional:
% - strWindType: Two testcases: step sweep 4 to 25 ms and normal dist. with
%   18 m/s mean (Default: Sweep)
% - loadData: load simulation output data if available instead of running
%   simulation (Default: 1)
% - figNo1: Number of figure (Default: 1)
% - useFASTForComparison: Uses FAST(1) or Simulink (0) PI data for comparison
%   (Default: 0, PI on the same simplified Simulink model as the MPC)
% - plantMdl: WECS model in the closed loop of both controllers: 'Simulink'
%   (simplified model Mdl_BianchiOL, Default), 'SimFilt' (simplified model
%   with the rotor-effective wind filter, Mdl_BianchiOL_WindFiltered, for
%   turbulent wind) or 'FAST' (FAST v8.16 S-function, input
%   files in FASTinputfiles; 'Sweep': FAST.fst, 18: FAST_NTW18.fst with a
%   TurbSim wind field, so the FAST 18 m/s wind differs from the
%   OutDataWind18NTW.mat wind, which is only used for initial values, the
%   stop time and the power reference) or 'SimulinkPaper' (as in the paper:
%   the test models as saved, i.e. the PI on Mdl_BianchiOL_WindFiltered
%   with the wind filter settings of the base workspace and the qLPV-MPC
%   on Mdl_BianchiOL)
% - ctrlOpts: qLPV-MPC variant, struct with the optional fields qMPC, rMPC,
%   pMPC (weights, replace the initWorkspace values during the run) and
%   name (suffix of the result files and figures). Default: the
%   initWorkspace weights, no suffix.
% - useFiltered: only 'SimulinkPaper': 1 (Default) uses the PI run with the
%   rotor-effective wind filter (OutTableTest2*_kVe0dot25_kVeTau0dot75.mat
%   of runCompareModels) if available, as in the paper; 2 (turbulent wind
%   only) additionally runs the qLPV-MPC on the same filtered plant
%   Mdl_BianchiOL_WindFiltered and schedules it with the filtered wind
%   (setFilteredWindQLMPC, result files *_windFilt); 0: no filter, both
%   controllers on the same plant Mdl_BianchiOL (the PI test model gets the
%   plant of the qLPV-MPC in memory, i.e. no filter and no projection factor
%   kWindEffFilt; PI result files *_noFilt)
% - strOut: suffix of the figure and LaTeX table files (e.g. '_old').
%   Default: no suffix.
%
% Evaluation window (figures, ratios, LaTeX table): stepped wind sweeps from
% 15 s, all other wind cases (turbulent 18 m/s) from 60 s, to the end of the
% wind data (sweep 400 s, turbulent 660 s); Simulink and FAST plant.

%% Handle optional inputs
% The default inputs are provided here.

if ~nargin || isempty(strWindType)
    strWindType = 'Sweep'; % 'NTW18';  
end

if nargin < 2 || isempty(loadData)
    loadData = 0; % load simulation output data if available;
end

if nargin < 3 || isempty(figNo1)
    figNo1 = 1;
end

if nargin < 4 || isempty(useFASTForComparison)
    useFASTForComparison = 0;
end

if nargin <5 || isempty(figDirStr)
    figDirStr = 'figDir7';
end

if nargin < 6 || isempty(useTitle)
    useTitle = [0,1]; %[0,0];
end

if nargin < 7 || isempty(plantMdl)
    plantMdl = 'Simulink';
end

if nargin < 8 || isempty(ctrlOpts)
    ctrlOpts = struct;
end

if nargin < 9 || isempty(useFiltered)
    useFiltered = 1;
end

if nargin < 10 || isempty(strOut)
    strOut = '';
end

%% Set path to directories

% Set path to inputdata directory
workDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(workDir);
dataInDir = fullfile( mainDir,'dataIn');
addpath(dataInDir);

% Set path to figure directory
figDir = fullfile(mainDir,figDirStr);
if ~isfolder(figDir)
    mkdir(figDir)
end

% Set path to data output directory (results of the paper models in dataOut/paper)
[~, ~, dataOutSub] = getWecsModelVersion;
dataDirOut = fullfile(workDir, 'dataOut', dataOutSub);
if ~isfolder(dataDirOut)
    mkdir(dataDirOut)
end

% Names of simulation output in correct order
varnames = {'Wind', 'RotSpeed', 'GenPwr', 'GenTq', 'BlPitch1', ...
    'NcIMUTAxs', 'NcIMUTAys'};

% qLPV-MPC weights: the controller reads qMPC, rMPC and pMPC from the base
% workspace; a variant in ctrlOpts replaces them during this run
strCtrl = '';
strTableFilt = ''; % suffix of the LaTeX table rows (filtered wind for the qLPV-MPC)
strPI = ''; % suffix of the PI result file (_noFilt: PI on Mdl_BianchiOL)
if isfield(ctrlOpts, 'name')
    strCtrl = ['_', ctrlOpts.name];
end
weightNames = {'qMPC', 'rMPC', 'pMPC'};
weightsOld = cell(size(weightNames));
for idxW = 1:numel(weightNames)
    weightsOld{idxW} = evalin('base', weightNames{idxW});
    if isfield(ctrlOpts, weightNames{idxW})
        assignin('base', weightNames{idxW}, ctrlOpts.(weightNames{idxW}));
    end
end
cleanupWeights = onCleanup(@() cellfun(@(n, v) assignin('base', n, v), weightNames, weightsOld));
q_ = evalin('base', 'qMPC');
r_ = evalin('base', 'rMPC');
p = evalin('base', 'pMPC')*q_;

% Controller description for the figure names: tracked outputs with weight
trackedNames = {'T_g', 'beta', 'omega_r'};
ctrlLabel = sprintf('qLMPC (tracking %s; q = [%s], r = [%s], p = %g*q)', ...
    strjoin(trackedNames(q_(1:3) > 0), ', '), num2str(q_), num2str(r_), p(1)/q_(1));

% Switch between Sweep and NTW18
if strcmp(strWindType,'Sweep') == 1 % sweep from 4 to 25 in steps
    outDataSimulationMat = 'OutDataSweep.mat'; 
    strFig = '';
elseif strcmp(strWindType,'Step') == 1 % sweep from 4 to 25 in steps
    [~, outDataSimulationMat] = getWecsModelVersion; % FAST step reference of the model version
    strFig = 'Step';
elseif isa(strWindType,'double')  % wind with average 18 m/s
    outDataSimulationMat = sprintf('OutDataWind%02dNTW.mat',strWindType);
    strFig = sprintf('NTW%02d',strWindType); %'NTW18';
else
    % EOG16mpers
    outDataSimulationMat = 'EOG16mpers.mat'; 
    strFig = 'EOG';
end





% Switch between the Simulink and the FAST WECS model
switch plantMdl
    case 'SimulinkPaper'
        simMdlname2 = 'test_SimulinkMdl2_Baseline';
        simMdlname = 'test_SimulinkMdl2_qLPVMPCbeta';
        FAST_InputFileName = '';
        strPlant = '';
        strFigPlant = 'Simulink';
        % WECS models as saved in the test models (discard edits in memory)
        wecsMdlCell = cell(1, 2);
        mdlCellPaper = {simMdlname2, simMdlname};
        for idxM = 1:2
            if bdIsLoaded(mdlCellPaper{idxM}), close_system(mdlCellPaper{idxM}, 0); end
            load_system(mdlCellPaper{idxM});
            blk = find_system(mdlCellPaper{idxM}, 'BlockType', 'ModelReference', 'Name', 'WECS Model');
            wecsMdlCell{idxM} = strrep(get_param(blk{1}, 'ModelName'), '.slx', '');
        end
        if useFiltered == 2 && isa(strWindType, 'double')
            % Turbulent wind: qLPV-MPC on the filtered plant of the PI,
            % scheduled with the filtered wind (model changed in memory only)
            oldBaseVars = setFilteredWindQLMPC(simMdlname);
            cleanupFilt = onCleanup(@() restoreBaseVars(oldBaseVars));
            cleanupMdlFilt = onCleanup(@() close_system(simMdlname, 0));
            wecsMdlCell{2} = 'Mdl_BianchiOL_WindFiltered, filtered wind';
            strCtrl = [strCtrl, '_windFilt'];
            strTableFilt = '_windFilt';
        end
        if useFiltered == 0
            % Both controllers on the same plant: the PI test model gets the
            % plant of the qLPV-MPC, Mdl_BianchiOL (no filter, no projection
            % factor kWindEffFilt; changed in memory only)
            blk = find_system(simMdlname2, 'BlockType', 'ModelReference', 'Name', 'WECS Model');
            set_param(blk{1}, 'ModelName', 'Mdl_BianchiOL');
            cleanupMdlPI = onCleanup(@() close_system(simMdlname2, 0));
            wecsMdlCell{1} = 'Mdl_BianchiOL';
            strPI = '_noFilt';
        end
        wecsLabel = sprintf('Mdl2 (PI: %s, qLMPC: %s)', wecsMdlCell{:});
    case {'Simulink', 'SimFilt'}
        simMdlname2 = 'test_SimulinkMdl2_Baseline';
        simMdlname = 'test_SimulinkMdl2_qLPVMPCbeta';
        FAST_InputFileName = '';
        if strcmp(plantMdl, 'Simulink')
            wecsMdl = 'Mdl_BianchiOL';
            strPlant = '';
            strFigPlant = 'Simulink';
        else
            % Rotor-effective wind in the plant aerodynamics, as in runCompareModels
            wecsMdl = 'Mdl_BianchiOL_WindFiltered';
            assignin('base', 'step_thresshold', 1);
            assignin('base', 'kVeFilt', 0.25);
            assignin('base', 'kVeFilt_tau', 0.75);
            strPlant = 'SimFilt';
            strFigPlant = 'SimulinkFiltered';
        end
        wecsLabel = ['Mdl2 (', wecsMdl, ')'];
        % Same plant for both controllers (changed in memory only)
        for mdlCell = {simMdlname2, simMdlname}
            load_system(mdlCell{1});
            blk = find_system(mdlCell{1}, 'BlockType', 'ModelReference', 'Name', 'WECS Model');
            set_param(blk{1}, 'ModelName', wecsMdl);
        end
        cleanupMdl = onCleanup(@() cellfun(@(m) close_system(m, 0), {simMdlname2, simMdlname}));
    case 'FAST'
        simMdlname2 = 'test_FAST_Baseline';
        simMdlname = 'test_FAST_qLPVMPCbeta.slx';
        if strcmp(strWindType,'Sweep') == 1
            FAST_InputFileName = fullfile(workDir,'FASTinputfiles','FAST.fst');
        elseif isequal(strWindType,18)
            FAST_InputFileName = fullfile(workDir,'FASTinputfiles','FAST_NTW18.fst');
        else
            error('No FAST input files for wind type %s.', string(strWindType));
        end
        strPlant = 'FAST';
        strFigPlant = 'FAST';
        [~, fstName] = fileparts(FAST_InputFileName);
        wecsLabel = ['FAST v8.16 (', fstName, '.fst)'];
    otherwise
        error('Unknown plant model %s (Simulink, SimulinkPaper, SimFilt or FAST).', plantMdl);
end

%% Load baseline and MPC data or run simulations

%For FASTtool simulation data
load(outDataSimulationMat ,'OutTable');

% Load data closed loop simulation PI
if useFASTForComparison % from FAST
    OutTableTest2 = OutTable;
else % from Simulink
    loadBaseline = 1; %always load if available
    matFileOutTableTest2 = fullfile(dataDirOut,['OutTableTest2',strPlant,strFig,strPI,'.mat']);
    % As in the paper: PI run with the rotor-effective wind filter
    % (runCompareModels, ..._kVe0dot25_kVeTau0dot75.mat) if available
    matFileOutTableTest2Filt = strrep(matFileOutTableTest2,'.mat','_kVe0dot25_kVeTau0dot75.mat');
    % (useFiltered = 2: the filter settings are set by setFilteredWindQLMPC,
    % so a missing file is simulated with the filter)
    if strcmp(plantMdl, 'SimulinkPaper') && ((useFiltered == 1 && exist(matFileOutTableTest2Filt,"file") == 2) || ...
            (useFiltered == 2 && isa(strWindType, 'double')))
        matFileOutTableTest2 = matFileOutTableTest2Filt;
    end
    OutTableTest2 = getSimulationOutputTable(matFileOutTableTest2,loadBaseline,OutTable,simMdlname2,[],FAST_InputFileName);
end

% Check that SI units are used (can be removed,
if mean(OutTableTest2.GenTq) <100 % protection against legacy data in kNm
    OutTableTest2.GenTq = OutTableTest2.GenTq*1000; % kNm -> Nm
end

if mean(OutTableTest2.RotSpeed) >1
    OutTableTest2.RotSpeed = OutTableTest2.RotSpeed *pi/30; %RPM -> rad/s
end

% Run or load qLPV controller
matFileOutTableTest1 = fullfile(dataDirOut,['OutTableMPC',strPlant,strFig,strCtrl,'.mat']);
[OutTableMPC, tictoc_LPVMPC,GenPwrRef] = getSimulationOutputTable(matFileOutTableTest1,loadData,OutTable,simMdlname,varnames,FAST_InputFileName);

%% Prepare information for plots and plot data

% Information for plot
figStr = ['Ref',strFigPlant,'ClosedLoop',strFig,'_beta',strCtrl,strOut];
strTableFilt = [strTableFilt, strOut]; % suffix of the LaTeX table file
% Figure name: controller and WECS model
if ischar(strWindType), strWind = strWindType; else, strWind = sprintf('%g m/s turbulent', strWindType); end
figLabel = sprintf('%s vs. PI | WECS: %s | wind: %s', ctrlLabel, wecsLabel, strWind);
meanOpt = strrep(sprintf(' mean cputime MPC: %1.2e s',mean(tictoc_LPVMPC.Data)),'e-0','e-');
% Evaluation window: stepped wind sweeps from 15 s (the start-up transient
% has decayed within the first 20 s wind level), otherwise from 60 s
if ischar(strWindType) && any(strcmp(strWindType, {'Sweep', 'Step'}))
    timeForPlot = [15, OutTable.Time(end)];
else
    timeForPlot = [60, OutTable.Time(end)];
end
DT = 0.008;
r = GenPwrRef;

% Plot and compare closed-loop results
plotComparison(OutTableTest2,OutTableMPC,q_,r_,p,DT,figStr,figDir,r,timeForPlot,meanOpt,figNo1,useTitle,strWindType,figLabel,strTableFilt)

function plotComparison(OutTableTest2,OutTableMPC,q_,r_,p,DT,figStr,figDir,r,timeForPlot,meanOpt,figNo1,useTitle,strWindType,figLabel,strTableFilt)
%plotComparison plots the time series obtained with baseline and qLPVMPC

%% Handle optional inputs
maxTime = min(height(OutTableTest2),height(OutTableMPC));
if nargin <10
    timeForPlot(1) = 0;
    timeForPlot(2) = maxTime;
end

%% Vector for plot: Time index for plot
timeVec = 0:DT: maxTime*DT - DT;

idxPlot = timeVec >= timeForPlot(1) & timeVec <= timeForPlot(2);
timeVecPlot = timeVec(idxPlot);
OutTableTest2Plot = OutTableTest2(idxPlot,:);
OutTableMPCPlot = OutTableMPC(idxPlot,:);

rData = squeeze(r.Data);
if numel(rData) == numel(timeVec)
    PGRef = rData(idxPlot);
else % reference logged at the controller sample time DTc
    PGRef = interp1(r.Time, rData, timeVecPlot(:), 'previous', 'extrap');
end
PGRef = PGRef(:);

%% Title strings
% Display optimization weights and calculate decrease in tower movement and
% power variance.

% String of optimization weights
qVec = deblank(sprintf('%d ',ceil(q_)));
rVec = deblank(sprintf('%d ', ceil(r_)));
pVec = deblank(sprintf('%d', p(1)/q_(1)));

% Improvement tower foreaft acceleration
varPI = var(OutTableTest2Plot.NcIMUTAxs);
varqLPV = var(OutTableMPCPlot.NcIMUTAxs);
titleStdPI = sprintf('%2.2e',varPI);
titleStdqLPV = sprintf('%2.2e',varqLPV);
ratioStd = sprintf('%2.1f ', varPI/varqLPV);

varPI_Pwr = var(PGRef - OutTableTest2Plot.GenPwr);
varqLPV_Pwr = var(PGRef - OutTableMPCPlot.GenPwr);


titleStdPI_Pwr = sprintf('%2.2e',varPI_Pwr);
titleStdqLPV_Pwr = sprintf('%2.2e',varqLPV_Pwr);
ratioStd_Pwr = sprintf('%2.1f ', varPI_Pwr/varqLPV_Pwr);

% --- Sample time ---
dt = mean(diff(timeVecPlot));

% --- Compute derivatives ---
dGenTq_PI   = diff(OutTableTest2Plot.GenTq/1000) / dt;  % kNm/s
dGenTq_qLPV = diff(OutTableMPCPlot.GenTq/1000)  / dt;
dPitch_PI   = diff(OutTableTest2Plot.BlPitch1)   / dt;  % deg/s
dPitch_qLPV = diff(OutTableMPCPlot.BlPitch1)     / dt;

% --- Define metrics ---
% Column 5: number format of the standard deviations in the LaTeX table
% (tower acceleration with four digits, all other signals with two digits)
metrics = {
    %  Label                      std_PI                                  std_qLPV                                unit             format
    'Tower fore-aft acc.',        std(OutTableTest2Plot.NcIMUTAxs),       std(OutTableMPCPlot.NcIMUTAxs),        'm/s$^2$',       '%.4f';
    'Power tracking error',       std(PGRef - OutTableTest2Plot.GenPwr),  std(PGRef - OutTableMPCPlot.GenPwr),   'kW',            '%.2f';
    'Generator torque rate',      std(dGenTq_PI),                         std(dGenTq_qLPV),                      'kNm/s',         '%.2f';
    'Blade pitch rate',           std(dPitch_PI),                         std(dPitch_qLPV),                      '\textdegree/s', '%.2f';
    };

% Treat numerically constant signals (e.g. constant torque above rated) as
% zero, so ratios become Inf/NaN instead of huge round-off values
metrics(:,2:3) = cellfun(@(x) x * (x > 1e-6), metrics(:,2:3), 'UniformOutput', false);

% --- Print std and variance ratios to console (for use in text) ---
fprintf('\n--- Std and variance ratios (PI/qLMPC, for text) ---\n');
for i = 1:size(metrics, 1)
    s_PI  = metrics{i,2};
    s_MPC = metrics{i,3};
    fmt   = metrics{i,5};
    fprintf(['%s: std_PI = ', fmt, ', std_qLMPC = ', fmt, ...
        ', std ratio = %.2f, var ratio = %.2f\n'], ...
        metrics{i,1}, s_PI, s_MPC, s_PI/s_MPC, (s_PI/s_MPC)^2);
end

% --- Write LaTeX table rows in thesis format: $std_PI / std_qLMPC = ratio$ ---
if isa(strWindType,'double')  % wind with average 18 m/s
    strWindType = sprintf('NTW%02d',strWindType); %'NTW18';
end

fileTable = ['results_table',strWindType,strTableFilt,'.tex']; % _windFilt: useFiltered = 2
fid = fopen(fileTable, 'w');
fprintf(fid, '%% Rows for table tab:ControllerComparisonMerged, wind case %s\n', strWindType);
fprintf(fid, '%% Standard deviations PI / qLMPC = ratio; ratios computed from unrounded values\n');
for i = 1:size(metrics, 1)
    s_PI  = metrics{i,2};
    s_MPC = metrics{i,3};
    fmt   = metrics{i,5};
    if isfinite(s_PI/s_MPC)
        cellStr = sprintf(['$', fmt, ' / ', fmt, ' = %.2f$'], s_PI, s_MPC, s_PI/s_MPC);
    else
        cellStr = sprintf(['$', fmt, ' / ', fmt, ' = {-}$'], s_PI, s_MPC);
    end
    fprintf(fid, '%s (%s) & %s \\\\\n', metrics{i,1}, metrics{i,4}, cellStr);
end
fclose(fid);
fprintf('Table rows written to %s\n', fileTable);
% axis for plot
yAxCell = {'Wind V (m/s)', 'Twr_{FA} (m/s^2)','GenPwr P_g (kW)'};

%% Plot
% Plot input wind, tower fore-aft acceleration and generator power

if useTitle(2) == 0 % either plot the
    figure(figNo1);

    % Plot wind
    ax1(1) = subplot(3,1,1);
    plot(timeVecPlot, OutTableMPCPlot.Wind);
    ylabel(yAxCell{1});
    axis tight; grid on;
    if useTitle(1)
        title({['Baseline (PI) vs. qLMPC; ', meanOpt],... %WT mdl: Model2
            ['q_ = [', qVec,'], r = [',rVec,'], p = ',pVec,'* q']})
    end

    % Tower fore-aft acceleration
    ax1(2) = subplot(3,1,2);
    plot(timeVecPlot, OutTableTest2Plot.NcIMUTAxs, timeVecPlot,OutTableMPCPlot.NcIMUTAxs,'-.');
    if useTitle(1)
        title([yAxCell{2}, ': var_{PI}: ',titleStdPI , ', var_{qLMPC}: ' ,titleStdqLPV,', ratio: ' ,ratioStd])
    end
    ylabel(yAxCell{2}); %'y_t [m/s^2]')
    axis tight; grid on;

    % Generator power
    ax1(3) = subplot(3,1,3);
    plot(timeVecPlot, OutTableTest2Plot.GenPwr/1000, timeVecPlot,OutTableMPCPlot.GenPwr/1000,'-.',timeVecPlot,PGRef/1000,'k--');
    if useTitle(1)
        title(['P_{g,ref}-P_g [kW]: var_{PI}: ',titleStdPI_Pwr, ', var_{qLMPC}: ' ,...
            titleStdqLPV_Pwr,', ratio: ' ,ratioStd_Pwr])
    end
    ylabel(strrep(yAxCell{3}, 'kW','MW'))
    axis tight; grid on;
    legend('PI','MPC','P_{g,ref}','Location','SouthEast')
    xlabel('Time (s)');

    % Link axes and set poisition
    linkaxes(ax1,'x')
    set(gcf,'Name', figLabel, 'NumberTitle', 'off')

    posDefault = [520   378   560   420]; %get(gcf, 'position');
    set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.1]);

    set(findall(gcf,'-property','FontSize'),'FontSize',12.5)
    set(findall(gcf,'-property','LineWidth'),'LineWidth',1.2)


    print(gcf,[fullfile(figDir,'aCmpCtrlSimulink_PI_MPC'),'_',figStr,'_',num2str(timeVecPlot(end))], '-dpng');
    print(gcf,[fullfile(figDir,'aCmpCtrlSimulink_PI_MPC'),'_',figStr,'_',num2str(timeVecPlot(end))], '-depsc');

    figNo1 = figNo1 + 1;

    if length(useTitle) == 2
        if (useTitle(2) == 0), return, end
    end
end

%% Figure for paper and dissertation
figure(figNo1+1)

% dc = [double('x'),776]; chX = char(dc);
%    yAxCell = {'Wind V (m/s)', 'GenTq T_g (kNm)', 'Pitch \beta (°)', 'RotSpd \omega_r (rpm)',...
%     'GenPwr P_g (MW)',['Twr_{FA}',chX,'_{t} (m/s^2)'], 'Twr_{SW} x_t (m/s^2)'};

% To do: Only names used for axis as x'' looks strange
% yAxCell = {'Wind V (m/s)', 'GenTq T_g (kNm)', 'Pitch \beta (°)', 'RotSpd \omega_r (rpm)',...
%  'TwrAcc_{FA} (m/s^2)', 'TwrAcc_{SW} (m/s^2)', 'GenPwr P_g (MW)'};

yAxCell = {'Wind V (m/s)', 'GenTq T_g (kNm)', 'Pitch \beta (°)', 'RotSpd \omega_r (rpm)',...
    'TwrAcc_{FA} (m/s^2)', 'TwrAcc_{SW} (m/s^2)', 'GenPwr P_g (MW)'};

yAxCell5 = {'Wind (m/s)', 'GenTq (kNm)', 'Pitch (°)', ...
    'TwrAcc_{FA} (m/s^2)', 'GenPwr (MW)',};


cl = lines;

titleStr = sprintf(['\\color[rgb]{%.3f %.3f %.3f}Baseline PI', ...
    '\\color{black} vs. ', ...
    '\\color[rgb]{%.3f %.3f %.3f}qLMPC', ...
    '\\color{black}; %s'], ...
    cl(1,1), cl(1,2), cl(1,3), ...
    cl(2,1), cl(2,2), cl(2,3), ...
    meanOpt);

%titleStr = ['Baseline PI vs. qLPV MPC; ', meanOpt];
nAx = 5;
tiledlayout(nAx,1,'TileSpacing','Compact');

if nAx ~= 7
    yAxCell = yAxCell5;
end

idxT = 1;

axPlotAll(idxT) = nexttile; %subplot(nAx,1,1);
plot(timeVecPlot,OutTableMPCPlot.Wind);
axis tight; grid on;
ylabel(yAxCell{idxT}); %'wind [m/s]')
title(titleStr);
idxT = idxT + 1;

axPlotAll(idxT) = nexttile; %axPlotAll(2) = subplot(nAx,1,2);
plot(timeVecPlot,OutTableTest2Plot.GenTq/1000,timeVecPlot,OutTableMPCPlot.GenTq/1000,'-.');
axis tight;
posAxis = axis;
axis([posAxis(1:2), min(43,posAxis(3)), 44]);
grid on; % axis tight;
ylabel(yAxCell{idxT}); %'GenTq T_g [Nm]')
idxT = idxT + 1;

axPlotAll(idxT) = nexttile; %axPlotAll(3) = subplot(nAx,1,3);
plot(timeVecPlot,OutTableTest2Plot.BlPitch1,timeVecPlot,OutTableMPCPlot.BlPitch1,'-.');
axis tight; grid on;
ylabel(yAxCell{idxT}); %'Pitch \beta [°]')
idxT = idxT + 1;

if nAx == 7
    axPlotAll(idxT) = nexttile; %axPlotAll(4) = subplot(nAx,1,4);
    plot(timeVecPlot,OutTableTest2Plot.RotSpeed * 30/pi,timeVecPlot,OutTableMPCPlot.RotSpeed,'-.');
    axis tight; grid on;
    ylabel(yAxCell{idxT}); %'RotSpd \omega_r [rpm]')
    idxT = idxT + 1;
end

axPlotAll(idxT) = nexttile; %axPlotAll(5) = subplot(nAx,1,5);
plot(timeVecPlot,OutTableTest2Plot.NcIMUTAxs,timeVecPlot,OutTableMPCPlot.NcIMUTAxs,'-.');
axis tight; grid on;
ylabel(yAxCell{idxT}); %'Twr_{FA} y_t[m/s^2]')
idxT = idxT + 1;

if nAx == 7
    axPlotAll(idxT) = nexttile; %axPlotAll(6) = subplot(nAx,1,6);
    plot(timeVecPlot,OutTableTest2Plot.NcIMUTAys,timeVecPlot,OutTableMPCPlot.NcIMUTAys,'-.');
    axis tight; grid on;
    ylabel(yAxCell{idxT}); %'Twr_{SW} x_t[m/s^2]')
    set(axPlotAll(6),'YLim', get(axPlotAll(5) ,'YLim'))
    idxT = idxT + 1;
end

axPlotAll(idxT) = nexttile; %axPlotAll(7) = subplot(nAx,1,7);
plot(timeVecPlot,OutTableTest2Plot.GenPwr/1000,timeVecPlot,OutTableMPCPlot.GenPwr/1000,'-.');
hold on; plot(timeVecPlot,PGRef/1000,'k--')
axis tight; grid on;
ylabel(yAxCell{idxT}); %'GenPwr P_g [MW]')

ax_pwr = gca;
text(ax_pwr, ...
    0.55, 0.17, 'P_{g,ref} (black dashed line)', ...
    ...  % 0.85, 0.17, 'P_{g,ref}',..
    'Units', 'normalized', ...
    'Interpreter', 'tex', ...
    'FontSize', 11, ...
    'BackgroundColor', [1 1 1 0.7], ...
    'Margin', 0.01, ...          % default is 3, reducing tightens the box
    'Color', 'k');

xlabel('Time (s)')
linkaxes(axPlotAll,'x');

set(gcf,'Name',figLabel,'NumberTitle','off')
posDefault = get(0,'DefaultFigurePosition');
set(gcf, 'position', [posDefault(1),posDefault(2) - posDefault(4)*0.7,posDefault(3),posDefault(4)*1.7]);

set(findall(gcf,'-property','FontSize'),'FontSize',11.5)
set(findall(gcf,'-property','LineWidth'),'LineWidth',0.75)

%set(groot,'defaultLineLineWidth',defaultLineWidth);

print(fullfile(figDir,['cmpCtrlTimeDomain_All_',figStr]), '-dpng');
print(fullfile(figDir,['cmpCtrlTimeDomain_All_',figStr]), '-depsc');

