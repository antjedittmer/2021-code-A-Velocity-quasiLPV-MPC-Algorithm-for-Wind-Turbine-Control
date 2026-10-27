%% Pareto comparison of qLMPC, IPC and AALC for the 18 m/s verification wind
% Evaluation criteria of the qLMPC, the IPC and the AALC (IPC + TEFC) for
% the turbulent 18 m/s wind of the second verification set (IEC class B,
% wind18_KlasseB_Verif2.bts = Eigene_Wind_Dateien/18mps.bts of the IPC/TEFC
% repository Combined-IPC-TEFC-for-Load-Reduction-of-Wind-Turbines),
% normalized with the PI baseline at the same pitch rate limit:
% - qLMPC on FAST v8.16 (designs flapNotch, flapNotchIPC, flapNotchSmooth12
%   and flapNotchSmooth12IPC of
%   qlmpcDesign.m, FASTinputfiles/V2B18/FAST_V2B18.fst) and the PI on FAST
% - qLMPC and PI on the simplified model Mdl2 (Mdl_BianchiOL, no model
%   mismatch; paper model version and qLMPC design 'paper' as in
%   runGeneratePicsPaper), with the hub height wind of a FAST run of the PI in
%   the same wind field (dataIn/OutDataWind18V2B.mat, generated here as
%   generateOutDataWind18NTW)
% - IPC (gains optimized with fmincon for each rate limit) and AALC (gains
%   optimized at 18 m/s for 8 deg/s) on OpenFAST v3.2.1, from the saved
%   results of checkIPC_WindSpeeds of the IPC/TEFC repository (folders
%   Main/Post_processing/IPCcheck_fixedRL<rate>_KlasseB_T60 and
%   IPCcheck_AALCfromIPC_RL08_KlasseB_T60, 18mps; runs with the FAST output
%   from 60 s)
% for the pitch rate limits 4, 8 and 13 deg/s (AALC: 8 deg/s only).
%
% Criteria (ratio to the PI, < 1 = reduction), evaluation from 60 s (the
% usual start-up time of the 660 s runs):
% - C_DEL_b: blade root out-of-plane DEL (RootMyc = RootMOoP, Woehler
%   exponent 10, sum of the three blades), C_DEL_t: tower base fore-aft
%   DEL (TwrBsMyt, exponent 4), C_DEL_tt: tower top fore-aft DEL
%   (YawBrMyp, exponent 4), all with calculateDamageBladesNew
% - C_P: mean of the generator power, C_dP: std of the generator power
% - C_Act: actuator power, mean over the blades and time of |pitch rate *
%   pitch moment RootMzb| (as makeAALC_SummaryTable)
% - proxies, also for Mdl2 (no blade and tower loads, no pitch moment):
%   C_TwrAcc (std of the tower top fore-aft acceleration NcIMUTAxs) and
%   C_PitchRate (mean of |pitch rate|)
% - C_ActEst: C_Act with the pitch moment of blade 1 estimated from the
%   pitch angle, its derivatives and the hub height wind (model fitted to
%   the FAST runs of the PI and flapNotch at 8 deg/s), for Mdl2 and, as a
%   check against C_Act, for the qLMPC on FAST
% Mdl2 gives C_P, C_dP, C_TwrAcc, C_PitchRate and C_ActEst only.
%
% The simulations run one after another (FAST about 4 min each) and are
% stored in NonLinMdl/dataOut/pareto (loaded if available, loadData = 1).
% The table is saved in figDirPareto/ParetoVerif2.mat and .csv; plot it
% with plotParetoVerif2. The data dictionaries are restored at the end.
clc; clear; close all;

Simulink.data.dictionary.closeAll('-discard');

%% Set path to initialization script initWorkspace and run it
addpath(genpath('NonLinMdl'));
addpath(genpath('LinMdl'));

% Pitch rate limits of the runs: data dictionaries restored at the end
copyDataDictionaries('backup');

initWorkspace;

%% Settings
loadData = 1; % load the stored simulation results if available
rateLimits = [4, 8, 13]; % deg/s
designNames = {'flapNotch', 'flapNotchIPC', 'flapNotchSmooth12', 'flapNotchSmooth12IPC'}; % qLMPC designs on FAST
tEvalStart = 60; % s

%% Directories
mainDir = fileparts(which('runParetoVerif2'));
workDir = fullfile(mainDir, 'NonLinMdl');
fstFile = fullfile(workDir, 'FASTinputfiles', 'V2B18', 'FAST_V2B18.fst');
windMat = fullfile(mainDir, 'dataIn', 'OutDataWind18V2B.mat');
dataDirOut = fullfile(workDir, 'dataOut', 'pareto');
if ~isfolder(dataDirOut), mkdir(dataDirOut); end
ipcDir = fullfile(fileparts(mainDir), 'Combined-IPC-TEFC-for-Load-Reduction-of-Wind-Turbines', ...
    'Main', 'Post_processing');
outDir = fullfile(mainDir, 'figDirPareto');
if ~isfolder(outDir), mkdir(outDir); end
cacheDir = fullfile(tempdir, 'qlmpc_pareto');
if ~isfolder(cacheDir), mkdir(cacheDir); end
Simulink.fileGenControl('set', 'CacheFolder', cacheDir, 'CodeGenFolder', cacheDir, 'createDir', true);

%% Wind of the verification set: hub height wind of a FAST run of the PI
% The FAST runs need the full field wind file
% NonLinMdl/FASTinputfiles/wind18_KlasseB_Verif2.bts (82 MB, not under
% version control). It is 18mps.bts of the IPC/TEFC repository
% https://github.com/antjedittmer/Combined-IPC-TEFC-for-Load-Reduction-of-Wind-Turbines
% (FAST_Sim_TEF/5MW-Turbine-Inputfiles/Sub-Inputfiles/WindData/Eigene_Wind_Dateien/18mps.bts);
% copy it there before new FAST runs.
if ~isfile(windMat)
    setPitchRateLimit('DD_CtrlBaseline.sldd', 'Control', 13);
    generateWindTable(fstFile, windMat, 660, DT);
end
load(windMat, 'OutTable');

%% Mdl2: PI and qLMPC on Mdl_BianchiOL
% As runGeneratePicsPaper: paper model version and qLMPC design 'paper'
% (first, the FAST runs below assign the qLMPC design parameters to the
% base workspace)
wecsModelVersion = 'paper';
designMPC = qlmpcDesign('paper');
fn = fieldnames(designMPC);
for idxF = 1:numel(fn)
    assignin('base', fn{idxF}, designMPC.(fn{idxF}));
end
ctrlDesign = 'paper';
updateDDMdl1(0.75, [], [], 0.9475);
mdlPI = 'test_SimulinkMdl2_Baseline';
mdlMPC = 'test_SimulinkMdl2_qLPVMPCbeta';
nR = numel(rateLimits);
runsMdl2 = cell(nR, 2);
for idxR = 1:nR
    rate = rateLimits(idxR);
    setPitchRateLimit('DD_CtrlBaseline.sldd', 'Control', rate);
    setPitchRateLimit('DD_Mdl1.sldd', 'controlValue', rate);
    mdlNames = {mdlPI, mdlMPC};
    for idxM = 1:2
        % Same plant for both controllers (in memory only)
        if bdIsLoaded(mdlNames{idxM}), close_system(mdlNames{idxM}, 0); end
        load_system(mdlNames{idxM});
        blk = find_system(mdlNames{idxM}, 'BlockType', 'ModelReference', 'Name', 'WECS Model');
        set_param(blk{1}, 'ModelName', 'Mdl_BianchiOL');
        strCtrl = {'PI', 'MPC'};
        matFile = fullfile(dataDirOut, sprintf('OutTableV2B18_rate%02d_%s.mat', rate, strCtrl{idxM}));
        runsMdl2{idxR, idxM} = compactMdl2Run(matFile, ...
            getSimulationOutputTable(matFile, loadData, OutTable, mdlNames{idxM}));
        close_system(mdlNames{idxM}, 0);
    end
end

%% FAST: PI and qLMPC designs
% As runGeneratePicsFASTthesis: FAST model version, data dictionaries as
% found (the paper version of the Mdl2 runs is undone)
copyDataDictionaries('reset');
wecsModelVersion = 'FAST';
ctrlNames = [{'PI'}, designNames];
runsFAST = cell(nR, numel(ctrlNames));
for idxR = 1:nR
    rate = rateLimits(idxR);
    setPitchRateLimit('DD_CtrlBaseline.sldd', 'Control', rate);
    setPitchRateLimit('DD_Mdl1.sldd', 'controlValue', rate);
    ipcMaxRateCtrl = rate; % IPC blade pitch (designs *IPC)
    for idxC = 1:numel(ctrlNames)
        matFile = fullfile(dataDirOut, sprintf('OutTable_V2B18_rate%02d_%s.mat', rate, ctrlNames{idxC}));
        runsFAST{idxR, idxC} = getFASTRun(matFile, loadData, ctrlNames{idxC}, OutTable, fstFile, rate);
    end
end

%% Evaluation criteria
% Model of the pitch moment RootMzb (blade 1) from the pitch angle, its
% derivatives and the hub height wind, fitted to the FAST runs of the PI and
% the qLMPC flapNotch at 8 deg/s: estimate of C_Act for Mdl2 (C_ActEst, also
% for the qLMPC on FAST as a check against C_Act)
windHub = OutTable.Wind1VelX;
idxR8 = find(rateLimits == 8);
fitRuns = [runsFAST(idxR8, 1), runsFAST(idxR8, strcmp(ctrlNames, 'flapNotch'))];
cMzb = fitPitchMoment(fitRuns, windHub, DT);
rows = {};
for idxR = 1:nR
    rate = rateLimits(idxR);
    % Mdl2 (proxies only)
    base = signalsMdl2(runsMdl2{idxR, 1}, DT);
    s = signalsMdl2(runsMdl2{idxR, 2}, DT);
    rows(end+1, :) = [criteriaRow("qLMPC", "Mdl2", rate, s, base, tEvalStart), ...
        {actEst(s, windHub, cMzb, DT, tEvalStart)/actEst(base, windHub, cMzb, DT, tEvalStart)}]; %#ok<SAGROW>
    % FAST v8.16
    base = signalsFAST(runsFAST{idxR, 1});
    for idxC = 2:numel(ctrlNames)
        s = signalsFAST(runsFAST{idxR, idxC});
        rows(end+1, :) = [criteriaRow("qLMPC " + ctrlNames{idxC}, "FAST", rate, s, base, tEvalStart), ...
            {actEst(s, windHub, cMzb, DT, tEvalStart)/actEst(base, windHub, cMzb, DT, tEvalStart)}]; %#ok<SAGROW>
    end
    % IPC (OpenFAST)
    S = load(fullfile(ipcDir, sprintf('IPCcheck_fixedRL%02d_KlasseB_T60', rate), '18mps', ...
        sprintf('TimeSeries_RL%02d.mat', rate)), 'TimeSeries');
    rows(end+1, :) = [criteriaRow("IPC", "OpenFAST", rate, signalsOpenFAST(S.TimeSeries.Optimum), ...
        signalsOpenFAST(S.TimeSeries.Baseline), tEvalStart), {NaN}]; %#ok<SAGROW>
end
S = load(fullfile(ipcDir, 'IPCcheck_AALCfromIPC_RL08_KlasseB_T60', '18mps', 'TimeSeries_RL08.mat'), 'TimeSeries');
rows(end+1, :) = [criteriaRow("AALC", "OpenFAST", 8, signalsOpenFAST(S.TimeSeries.Optimum), ...
    signalsOpenFAST(S.TimeSeries.Baseline), tEvalStart), {NaN}];
ParetoTable = cell2table(rows, 'VariableNames', {'Controller', 'Plant', 'RateLimit_degps', ...
    'C_DEL_b', 'C_DEL_t', 'C_DEL_tt', 'C_P', 'C_dP', 'C_Act', 'C_TwrAcc', 'C_PitchRate', 'C_ActEst'});
ParetoTable = sortrows(ParetoTable, {'Controller', 'RateLimit_degps'});
disp(ParetoTable);
save(fullfile(outDir, 'ParetoVerif2.mat'), 'ParetoTable', 'tEvalStart');
writetable(ParetoTable, fullfile(outDir, 'ParetoVerif2.csv'));
fprintf('Saved %s\n', fullfile(outDir, 'ParetoVerif2.mat'));

%% Restore the data dictionaries
copyDataDictionaries('restore');

%% Local functions
function copyDataDictionaries(direction)
% copyDataDictionaries copies the data dictionaries changed by the runs to
% tempdir ('backup') or back ('restore', as runGeneratePicsPaper; 'reset':
% back, but the copy is kept)
ddFiles = cellfun(@which, {'DD_test.sldd', 'DD_CtrlBaseline.sldd', 'DD_Mdl1.sldd'}, 'UniformOutput', false);
for idxD = 1:numel(ddFiles)
    [~, name, ext] = fileparts(ddFiles{idxD});
    ddBackup = fullfile(tempdir, ['runParetoVerif2_', name, ext, '.bak']);
    if strcmp(direction, 'backup')
        copyfile(ddFiles{idxD}, ddBackup, 'f');
    else
        if idxD == 1
            bdclose('all');
            Simulink.data.dictionary.closeAll('-discard');
        end
        copyfile(ddBackup, ddFiles{idxD}, 'f');
        if strcmp(direction, 'restore')
            delete(ddBackup);
        end
    end
end
end

function setPitchRateLimit(DDName, entryName, maxRate)
% setPitchRateLimit sets the pitch rate limits +/-maxRate (deg/s) of the
% entry entryName in the data dictionary DDName (saved only if changed)
dictObj = Simulink.data.dictionary.open(DDName);
entryObj = getEntry(getSection(dictObj, 'Design Data'), entryName);
val = getValue(entryObj);
if val.Pitch.Maxrate ~= maxRate || val.Pitch.Minrate ~= -maxRate
    val.Pitch.Maxrate = maxRate;
    val.Pitch.Minrate = -maxRate;
    setValue(entryObj, val);
    saveChanges(dictObj);
    fprintf('%s: pitch rate limit set to %g deg/s\n', DDName, maxRate);
end
end

function generateWindTable(fstFile, outMat, tEnd, DT)
% generateWindTable runs the PI on FAST in the wind field of fstFile and
% saves the FAST output channels of OutList_FAST.txt as OutTable (as
% generateOutDataWind18NTW: the run is repeated with the hub height wind of
% the first run, which sets the initial values of the test model)
inDir = fileparts(fstFile);
simMdlname = 'test_FAST_Baseline';
channels = cellstr(readlines(fullfile(inDir, 'OutList_FAST.txt'), 'EmptyLineRule', 'skip'));
channels = channels(:)';
time = (0:DT:tEnd)';
windSeed = 18*ones(size(time)); % first run: mean wind speed
for idxRun = 1:2
    OutTable = table(time, windSeed, 'VariableNames', {'Time', 'Wind1VelX'});
    assignin('base', 'OutTable', OutTable);
    assignin('base', 'FAST_InputFileName', fstFile);
    assignin('base', 'TMax', tEnd);
    assignin('base', 'OutList', channels');
    load_system(simMdlname);
    fprintf('FAST run %d of 2 (%s, %g s) ...\n', idxRun, fstFile, tEnd);
    sim(simMdlname, 'SrcWorkspace', 'base');
    close_system(simMdlname, 0);
    [~, fstName] = fileparts(fstFile);
    OutTableFAST = readFASTOut(fullfile(inDir, [fstName, '.SFunc.out']), channels);
    assert(abs(OutTableFAST.Time(end) - tEnd) < DT/2, 'FAST run ended at %g s', OutTableFAST.Time(end));
    windSeed = interp1(OutTableFAST.Time, OutTableFAST.Wind1VelX, time, 'linear', 'extrap');
end
OutTable = OutTableFAST;
OutTable.Time = time; % exact time grid
save(outMat, 'OutTable');
fprintf('Saved %s: wind mean %.2f m/s, std %.2f m/s\n', outMat, mean(OutTable.Wind1VelX), std(OutTable.Wind1VelX));
end

function R = getFASTRun(matFile, loadData, ctrlName, OutTable, fstFile, maxRate)
% getFASTRun loads the stored closed-loop FAST run or runs the simulation
% (as runGeneratePicsFASTthesis, with the FAST output channels of the
% evaluation criteria; stored runs without the tower top moment are repeated)
if loadData && isfile(matFile)
    tmp = load(matFile, 'R');
    if ~isempty(tmp.R.FASTOut) && ismember('YawBrMyp', tmp.R.FASTOut.Properties.VariableNames)
        R = tmp.R;
        fprintf('Loaded %s\n', matFile);
        return
    end
end

% Controller: PI baseline or qLPV-MPC design (parameters in the base workspace)
if strcmp(ctrlName, 'PI')
    design = [];
    simMdlname = 'test_FAST_Baseline';
else
    design = qlmpcDesign(ctrlName);
    fn = fieldnames(design);
    for k = 1:numel(fn)
        assignin('base', fn{k}, design.(fn{k}));
    end
    assignin('base', 'ctrlDesign', ctrlName);
    simMdlname = ['test_FAST_qLPVMPC', design.ctrlVariant];
end
assignin('base', 'ipcMaxRateCtrl', maxRate);

% FAST input file copy with its own root name (outputs FAST_exp_pareto.*)
fstRun = fullfile(fileparts(fstFile), 'FAST_exp_pareto.fst');
copyfile(fstFile, fstRun);
assignin('base', 'OutTable', OutTable);
assignin('base', 'FAST_InputFileName', fstRun);
assignin('base', 'TMax', OutTable.Time(end));
assignin('base', 'OutList', cellstr(readlines(fullfile(fileparts(fstFile), 'OutList_FAST.txt'), 'EmptyLineRule', 'skip')));

fprintf('Running %s with %s at %g deg/s (%s) ...\n', simMdlname, ctrlName, maxRate, fstFile);
R = struct('ctrlName', ctrlName, 'design', design, 'status', "ok", 'msg', "", ...
    'maxRate', maxRate, 'fstFile', fstFile, 'FASTOut', table());
t0 = tic;
try
    load_system(simMdlname);
    out = sim(simMdlname, 'SrcWorkspace', 'base', 'CaptureErrors', 'on');
    if ~isempty(out.ErrorMessage)
        R.status = "failed";
        R.msg = string(out.ErrorMessage);
    end
catch ME
    R.status = "failed";
    R.msg = string(ME.message);
end
R.wallTime_s = toc(t0);
close_system(simMdlname, 0);
outFile = fullfile(fileparts(fstRun), 'FAST_exp_pareto.SFunc.out');
if isfile(outFile)
    R.FASTOut = readFASTOut(outFile, {'Time', 'RootMOoP1', 'RootMOoP2', 'RootMOoP3', ...
        'BlPitch1', 'BlPitch2', 'BlPitch3', 'RootMzb1', 'RootMzb2', 'RootMzb3', ...
        'TwrBsMyt', 'YawBrMyp', 'NcIMUTAxs', 'GenPwr'});
end
delete(fullfile(fileparts(fstRun), 'FAST_exp_pareto.*'));
if R.status ~= "ok"
    warning('%s: %s', ctrlName, R.msg);
end
save(matFile, 'R');
fprintf('Saved %s (%.0f s)\n', matFile, R.wallTime_s);
end

function T = readFASTOut(outFile, channels)
% readFASTOut reads the channels of a FAST text output file (tab delimited,
% header with the channel names followed by a line with the units)
lines = readlines(outFile);
idxHead = find(startsWith(strtrim(lines), 'Time'), 1);
names = strtrim(split(lines(idxHead), char(9)));
dataLines = lines(idxHead + 2:end);
dataLines = dataLines(contains(dataLines, char(9))); % without empty lines at the end
data = str2double(split(strtrim(dataLines), char(9)));
data = data(all(isfinite(data), 2), :);
[isOut, idxCol] = ismember(channels, names);
assert(all(isOut), 'Channels missing in %s: %s', outFile, strjoin(channels(~isOut), ', '));
T = array2table(data(:, idxCol), 'VariableNames', channels);
end

function T = compactMdl2Run(matFile, T)
% compactMdl2Run keeps only the columns of a Mdl2 run used by signalsMdl2
% and stores the reduced table in matFile (without the qLPV-MPC timing)
keep = {'BlPitch1', 'GenPwr', 'NcIMUTAxs'};
if width(T) > numel(keep)
    T = T(:, keep);
    [~, strTable] = fileparts(matFile);
    S.(strTable) = T;
    save(matFile, '-struct', 'S');
    fprintf('Reduced %s to %s\n', matFile, strjoin(keep, ', '));
end
end

function c = fitPitchMoment(runs, windHub, DT)
% fitPitchMoment least squares fit of the pitch moment RootMzb1 of the FAST
% runs (struct cell array) to the regressors of pitchMomentRegressors
X = [];
y = [];
for k = 1:numel(runs)
    F = runs{k}.FASTOut;
    X = [X; pitchMomentRegressors(F.BlPitch1, windHub, DT)]; %#ok<AGROW>
    y = [y; F.RootMzb1]; %#ok<AGROW>
end
c = X\y;
fprintf('Pitch moment model: R^2 = %.2f\n', 1 - var(y - X*c)/var(y));
end

function X = pitchMomentRegressors(pitch, windHub, DT)
% pitchMomentRegressors quadratic in pitch angle and wind, pitch rate and
% pitch acceleration (pitch in deg, time step DT)
X = [ones(size(pitch)), pitch, pitch.^2, windHub, windHub.^2, pitch.*windHub, ...
    [0; diff(pitch)]/DT, [0; 0; diff(pitch, 2)]/DT^2];
end

function a = actEst(s, windHub, c, DT, tStart)
% actEst actuator power with the estimated pitch moment, blade 1
pitch = s.pitch(:, 1);
mzb = pitchMomentRegressors(pitch, windHub, DT)*c;
idx = find(s.t >= tStart);
a = mean(abs(diff(pitch(idx))/DT.*mzb(idx(1:end-1))));
end

function s = signalsMdl2(T, DT)
% Signals of a Mdl2 run (collective pitch, no loads)
s.t = (0:height(T) - 1)'*DT;
s.pitch = T.BlPitch1;
s.genPwr = T.GenPwr;
s.twrAcc = T.NcIMUTAxs;
end

function s = signalsFAST(R)
% Signals of a FAST v8.16 run (FAST output file)
F = R.FASTOut;
assert(R.status == "ok", '%s: %s', R.ctrlName, R.msg);
s.t = F.Time;
s.pitch = [F.BlPitch1, F.BlPitch2, F.BlPitch3];
s.genPwr = F.GenPwr;
s.twrAcc = F.NcIMUTAxs;
s.myc = [F.RootMOoP1, F.RootMOoP2, F.RootMOoP3]; % RootMOoP = RootMyc
s.mzb = [F.RootMzb1, F.RootMzb2, F.RootMzb3];
s.twrMy = F.TwrBsMyt;
s.twrTopMy = F.YawBrMyp;
end

function s = signalsOpenFAST(ts)
% Signals of an OpenFAST run of the IPC/TEFC repository (TimeSeries struct)
OutList = strtrim(string(ts.OutList));
ch = @(name) ts.OutData(:, OutList == name);
s.t = ch("Time");
s.pitch = [ch("BldPitch1"), ch("BldPitch2"), ch("BldPitch3")];
s.genPwr = ch("GenPwr");
s.twrAcc = ch("NcIMUTAxs");
s.myc = [ch("RootMyc1"), ch("RootMyc2"), ch("RootMyc3")];
s.mzb = [ch("RootMzb1"), ch("RootMzb2"), ch("RootMzb3")];
s.twrMy = ch("TwrBsMyt");
s.twrTopMy = ch("YawBrMyp");
end

function row = criteriaRow(ctrl, plant, rate, s, b, tStart)
% criteriaRow evaluation criteria of the run s, normalized with the PI run b
c = criteria(s, tStart);
cB = criteria(b, tStart);
row = {ctrl, plant, rate, c.delB/cB.delB, c.delT/cB.delT, c.delTT/cB.delTT, c.meanP/cB.meanP, c.stdP/cB.stdP, ...
    c.act/cB.act, c.stdTwrAcc/cB.stdTwrAcc, c.pitchRate/cB.pitchRate};
end

function c = criteria(s, tStart)
% criteria in the evaluation window t >= tStart (loads NaN without load signals)
idx = s.t >= tStart;
t = s.t(idx);
dt = diff(t);
rate = diff(s.pitch(idx, :))./dt;
c.meanP = mean(s.genPwr(idx));
c.stdP = std(s.genPwr(idx));
c.stdTwrAcc = std(s.twrAcc(idx));
c.pitchRate = mean(abs(rate), 'all');
[c.delB, c.delT, c.delTT, c.act] = deal(NaN);
if isfield(s, 'myc')
    c.delB = calculateDamageBladesNew(s.myc(idx, :), 10); % sum of the blades
    c.delT = calculateDamageBladesNew(s.twrMy(idx), 4);
    c.delTT = calculateDamageBladesNew(s.twrTopMy(idx), 4);
    mzb = s.mzb(idx, :);
    c.act = mean(abs(rate.*mzb(1:end-1, :)), 'all');
end
end
