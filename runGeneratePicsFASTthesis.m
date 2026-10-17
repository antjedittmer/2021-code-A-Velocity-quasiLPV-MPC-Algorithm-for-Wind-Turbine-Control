%% Script to compare the PI and the most promising qLPV-MPC designs on FAST
% Closed-loop runs of the baseline PI controller and the stored qLPV-MPC
% designs (subfunctionsMPC/qlmpcDesign.m) with FAST v8.16 as the WECS model
% for the turbulent 18 m/s wind (FAST_NTM18B.fst: TurbSim, IEC class B, NTM,
% shear 0.2, 660 s, the wind field of dataIn/OutDataWind18NTW.mat) and the wind sweep
% (FAST.fst). Pitch rate limit of both controllers: 13 deg/s.
%
% The simulation results are stored in NonLinMdl/dataOut/FASTthesis and
% loaded if available (loadData = 1); set loadData = 0 to rerun FAST.
% FAST must not run in a second MATLAB session at the same time (shared
% output files and data dictionaries).
%
% Outputs:
% - time series figures (PI and qLMPC designs) per wind case in figDirFASTthesis
% - comparison table per wind case: console, figDirFASTthesis/
%   resultsFASTthesis_<wind>.csv and .tex, NonLinMdl/dataOut/FASTthesis/
%   resultsFASTthesis.mat (variable resultsFAST)
%
% Metrics (evaluation window from 60 s: 60 to 660 s for the 18 m/s wind, 60 s to
% the end of the data for the sweep): std
% of the rotor speed, the power tracking error |P_g,ref - P_g|, the tower
% fore-aft acceleration, the generator torque rate, the (collective) pitch
% rate and the blade 1 pitch rate; from the FAST output file the std of the
% blade root out-of-plane moments (mean of the three blades), the blade root
% DEL (Woehler exponent 10, maximum of the three blades) and the tower base
% fore-aft DEL (TwrBsMyt, exponent 4), as in thesis section 4.2. The ratio
% PI/qLMPC is above 1 if the qLMPC design is better.
clc; clear; close all;

Simulink.data.dictionary.closeAll('-discard');

%% Set path to initialization script initWorkspace and run it
addpath(genpath('NonLinMdl'));
addpath(genpath('LinMdl'));

initWorkspace;

%% Settings
loadData = 1; % load the stored simulation results if available
maxRate = 13; % pitch rate limit (deg/s) for PI, qLMPC and the IPC blade pitch
ipcMaxRateCtrl = maxRate;

% qLPV-MPC designs (qlmpcDesign.m, results in NonLinMdl/experiments/EXPERIMENTS.md)
% - 'beta': original design (paper/thesis), tracks Tg and the pitch
% - 'omega': rotor speed tracking, lowest pitch activity and tower acceleration
% - 'flapNotch': tower/blade flap model with notch filters, best speed and power
% - 'flapNotchSmooth12': DTc = 0.12 s, smooth pitch reference, balanced
% - 'flapNotchSmooth12Tower': as above with the tower velocity weight
% - 'flapNotchIPC': flapNotch with the DLR individual pitch controller
%   (ipc_dlr, adapted to FAST: gain set 2 scaled by 0.3, no tower term)
designNames = {'beta', 'omega', 'flapNotch', 'flapNotchSmooth12', 'flapNotchSmooth12Tower', ...
    'flapNotchIPC'};

% Wind cases: name, reference data (initial values, stop time), FAST input file,
% label. NTM18B: the same class B wind field as OutDataWind18NTW.mat (660 s);
% the earlier case NTW18 (FAST_NTW18.fst, class A, 600 s) is no longer run.
windCases = {'NTM18B', 'OutDataWind18NTW.mat', 'FAST_NTM18B.fst', '18 m/s turbulent (class B)';
    'Sweep', 'OutDataSweep.mat', 'FAST.fst', 'sweep'};

tEvalStart = 60; % evaluation from 60 s to the end of the wind data (s)
figDirStr = 'figDirFASTthesis';

%% Directories
mainDir = fileparts(mfilename('fullpath'));
workDir = fullfile(mainDir, 'NonLinMdl');
dataDirOut = fullfile(workDir, 'dataOut', 'FASTthesis');
if ~isfolder(dataDirOut), mkdir(dataDirOut); end
figDir = fullfile(mainDir, figDirStr);
if ~isfolder(figDir), mkdir(figDir); end
cacheDir = fullfile(tempdir, 'qlmpc_FASTthesis');
if ~isfolder(cacheDir), mkdir(cacheDir); end
Simulink.fileGenControl('set', 'CacheFolder', cacheDir, 'CodeGenFolder', cacheDir, 'createDir', true);

%% Pitch rate limit: qLMPC (DD_Mdl1.sldd) and PI (DD_CtrlBaseline.sldd)
setPitchRateLimit('DD_Mdl1.sldd', 'controlValue', maxRate);
setPitchRateLimit('DD_CtrlBaseline.sldd', 'Control', maxRate);

%% Closed-loop runs and comparison
ctrlNames = [{'PI'}, designNames];
resultsFAST = struct;
figNo = 300;
for idxW = 1:size(windCases, 1)
    windName = windCases{idxW, 1};
    S = load(fullfile(mainDir, 'dataIn', windCases{idxW, 2}), 'OutTable');
    OutTable = S.OutTable;
    fstFile = fullfile(workDir, 'FASTinputfiles', windCases{idxW, 3});
    timeEval = [tEvalStart, OutTable.Time(end)]; % evaluation window (s)
    timeEvalCase.(windName) = timeEval;

    % Run or load the simulations
    runs = cell(size(ctrlNames));
    for idxC = 1:numel(ctrlNames)
        matFile = fullfile(dataDirOut, sprintf('OutTable_%s_%s.mat', windName, ctrlNames{idxC}));
        runs{idxC} = getFASTRun(matFile, loadData, ctrlNames{idxC}, OutTable, fstFile, maxRate);
    end

    % Power reference of the qLMPC (depends on the measured wind only, the
    % FAST wind field is the same for all controllers); used for all runs
    % (taken from the longest run, an aborted run has a shorter reference)
    refEnd = -Inf(size(runs));
    for idxC = 1:numel(runs)
        if ~isempty(runs{idxC}.GenPwrRef)
            refEnd(idxC) = runs{idxC}.GenPwrRef(end, 1);
        end
    end
    if all(isinf(refEnd))
        error('No power reference logged for wind case %s.', windName);
    end
    [~, idxRef] = max(refEnd);
    GenPwrRef = runs{idxRef}.GenPwrRef;

    % Metrics
    resultsFAST.(windName) = evaluateRuns(runs, ctrlNames, GenPwrRef, timeEval);
    printTable(resultsFAST.(windName), windCases{idxW, 4}, maxRate);
    writeTable(resultsFAST.(windName), figDir, windName, windCases{idxW, 4}, maxRate);

    % Time series
    plotRuns(runs, ctrlNames, GenPwrRef, timeEval, figNo, figDir, windName, windCases{idxW, 4}, fstFile, maxRate);
    figNo = figNo + 1;
end
save(fullfile(dataDirOut, 'resultsFASTthesis.mat'), 'resultsFAST', 'designNames', 'windCases', 'maxRate', 'timeEvalCase');

%% Local functions
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

function R = getFASTRun(matFile, loadData, ctrlName, OutTable, fstFile, maxRate)
% getFASTRun loads the stored closed-loop FAST run or runs the simulation.
% R: struct with the fields OutTableRun (table at DT), GenPwrRef ([time,
% value]), tictoc (qLMPC computation time per step), design (qlmpcDesign
% struct, empty for the PI), status, msg, maxRate, wallTime_s, FASTOut (table
% of FAST output channels: blade root moments, blade 1 pitch, azimuth,
% tower base moment)
if loadData && isfile(matFile)
    tmp = load(matFile, 'R');
    % runs without any data or without the FAST output channels are repeated
    if ~isempty(tmp.R.OutTableRun) && isfield(tmp.R, 'FASTOut') && ~isempty(tmp.R.FASTOut)
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

% FAST input file copy with its own root name (outputs FAST_exp_thesis.*)
fstRun = fullfile(fileparts(fstFile), 'FAST_exp_thesis.fst');
copyfile(fstFile, fstRun);
assignin('base', 'OutTable', OutTable);
assignin('base', 'FAST_InputFileName', fstRun);
assignin('base', 'TMax', OutTable.Time(end));
assignin('base', 'OutList', cellstr(readlines(fullfile(fileparts(fstFile), 'OutList_FAST.txt'), 'EmptyLineRule', 'skip')));

fprintf('Running %s with %s (%s) ...\n', simMdlname, ctrlName, fstFile);
R = struct('ctrlName', ctrlName, 'design', design, 'status', "ok", 'msg', "", ...
    'maxRate', maxRate, 'fstFile', fstFile, 'OutTableRun', table(), 'GenPwrRef', [], 'tictoc', [], ...
    'FASTOut', table());
t0 = tic;
try
    load_system(simMdlname);
    % CaptureErrors: keep the logged data up to an abort (e.g. FAST divergence)
    out = sim(simMdlname, 'SrcWorkspace', 'base', 'CaptureErrors', 'on');
    if ~isempty(out.ErrorMessage)
        R.status = "failed";
        R.msg = string(out.ErrorMessage);
    end
    % Logged signals (the PI model logs OutDataTest only)
    logged = out.who;
    varnames = {'Wind', 'RotSpeed', 'GenPwr', 'GenTq', 'BlPitch1', 'NcIMUTAxs', 'NcIMUTAys'};
    if any(strcmp(logged, 'OutDataTest'))
        D = out.OutDataTest;
        R.OutTableRun = array2table(D, 'VariableNames', varnames(1:size(D, 2)));
    end
    if any(strcmp(logged, 'GenPwrRef'))
        R.GenPwrRef = [out.GenPwrRef.Time(:), squeeze(out.GenPwrRef.Data(:))];
    end
    if any(strcmp(logged, 'tictoc_LPVMPC'))
        R.tictoc = squeeze(out.tictoc_LPVMPC.Data);
    end
catch ME
    R.status = "failed";
    R.msg = string(ME.message);
end
R.wallTime_s = toc(t0);
close_system(simMdlname, 0);
outFile = fullfile(fileparts(fstRun), 'FAST_exp_thesis.SFunc.out');
if isfile(outFile)
    R.FASTOut = readFASTOut(outFile, {'Time', 'RootMOoP1', 'RootMOoP2', 'RootMOoP3', ...
        'BlPitch1', 'Azimuth', 'TwrBsMyt'});
end
delete(fullfile(fileparts(fstRun), 'FAST_exp_thesis.*'));
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
[~, idxCol] = ismember(channels, names);
T = array2table(data(:, idxCol), 'VariableNames', channels);
end

function D = damageEquivalentLoad(x, t, m)
% damageEquivalentLoad DEL of the signal x for the Woehler exponent m with
% rainflow counting, equivalent number of cycles = duration in s (1 Hz)
c = rainflow(x, t);
D = (sum(c(:, 1).*c(:, 2).^m)/(t(end) - t(1)))^(1/m);
end

function [T, t] = runTable(R, DT)
% runTable returns the run's table in SI-like units (GenTq in Nm, RotSpeed
% in rpm) and its time vector
T = R.OutTableRun;
if ~isempty(T) && mean(T.GenTq) < 100 % legacy data in kNm
    T.GenTq = T.GenTq*1000;
end
t = (0:height(T) - 1)'*DT;
end

function Tres = evaluateRuns(runs, ctrlNames, GenPwrRef, timeEval)
% evaluateRuns computes the metrics of all runs in the evaluation window
DT = 0.008;
nC = numel(runs);
[RotSpdStd, RotSpdMax, PwrErrStd, TwrFAStd, GenTqRateStd, PitchRateStd, ...
    BlPitchRateStd, BlMOoPStd, BlDEL, TwrDEL, tCpuMean, simulated] = deal(nan(nC, 1));
status = strings(nC, 1);
for idxC = 1:nC
    [T, t] = runTable(runs{idxC}, DT);
    status(idxC) = runs{idxC}.status;
    simulated(idxC) = height(T)*DT;
    if isempty(T), continue, end
    idx = t >= timeEval(1) & t <= timeEval(2);
    PGRef = interp1(GenPwrRef(:, 1), GenPwrRef(:, 2), t(idx), 'previous', 'extrap');
    RotSpdStd(idxC) = std(T.RotSpeed(idx));
    RotSpdMax(idxC) = max(T.RotSpeed(idx));
    PwrErrStd(idxC) = std(abs(PGRef - T.GenPwr(idx)));
    TwrFAStd(idxC) = std(T.NcIMUTAxs(idx));
    GenTqRateStd(idxC) = std(diff(T.GenTq(idx)/1000)/DT);
    PitchRateStd(idxC) = std(diff(T.BlPitch1(idx))/DT);
    if ~isempty(runs{idxC}.tictoc)
        tCpuMean(idxC) = mean(runs{idxC}.tictoc);
    end
    % Blade and tower loads from the FAST output file
    F = runs{idxC}.FASTOut;
    if ~isempty(F)
        idxF = F.Time >= timeEval(1) & F.Time <= timeEval(2);
        tF = F.Time(idxF);
        BlPitchRateStd(idxC) = std(diff(F.BlPitch1(idxF))./diff(tF));
        Myb = [F.RootMOoP1(idxF), F.RootMOoP2(idxF), F.RootMOoP3(idxF)];
        BlMOoPStd(idxC) = mean(std(Myb));
        BlDEL(idxC) = max([damageEquivalentLoad(Myb(:, 1), tF, 10), ...
            damageEquivalentLoad(Myb(:, 2), tF, 10), damageEquivalentLoad(Myb(:, 3), tF, 10)]);
        TwrDEL(idxC) = damageEquivalentLoad(F.TwrBsMyt(idxF), tF, 4);
    end
end
% Numerically constant signals (constant torque above rated) count as zero
GenTqRateStd(GenTqRateStd < 1e-6) = 0;
Tres = table(RotSpdStd, RotSpdMax, PwrErrStd, TwrFAStd, GenTqRateStd, PitchRateStd, ...
    BlPitchRateStd, BlMOoPStd, BlDEL, TwrDEL, tCpuMean, simulated, status, 'RowNames', ctrlNames(:));
% Ratios PI/qLMPC (above 1: qLMPC better)
metricNames = {'RotSpdStd', 'PwrErrStd', 'TwrFAStd', 'GenTqRateStd', 'PitchRateStd', ...
    'BlPitchRateStd', 'BlMOoPStd', 'BlDEL', 'TwrDEL'};
for k = 1:numel(metricNames)
    Tres.(['ratio_', metricNames{k}]) = Tres.(metricNames{k})(1)./Tres.(metricNames{k});
end
end

function printTable(Tres, windLabel, maxRate)
% printTable prints the comparison table to the console
fprintf('\n=== FAST v8.16, %s, pitch rate limit %g deg/s ===\n', windLabel, maxRate);
fprintf('%-24s %8s %8s %8s %8s %9s %9s | %-34s\n', 'Controller', 'RotSpd', 'RotMax', 'PwrErr', 'TwrFA', ...
    'TqRate', 'PtchRate', 'std ratio PI/qLMPC: spd pwr twr pitch');
fprintf('%-24s %8s %8s %8s %8s %9s %9s |\n', '', '(rpm)', '(rpm)', '(kW)', '(m/s^2)', '(kNm/s)', '(deg/s)');
for k = 1:height(Tres)
    fprintf('%-24s %8.3f %8.2f %8.1f %8.3f %9.2f %9.2f | %5.2f %5.2f %5.2f %5.2f  %s\n', Tres.Row{k}, ...
        Tres.RotSpdStd(k), Tres.RotSpdMax(k), Tres.PwrErrStd(k), Tres.TwrFAStd(k), Tres.GenTqRateStd(k), ...
        Tres.PitchRateStd(k), Tres.ratio_RotSpdStd(k), Tres.ratio_PwrErrStd(k), Tres.ratio_TwrFAStd(k), ...
        Tres.ratio_PitchRateStd(k), Tres.status(k));
end
fprintf('\n%-24s %9s %9s %9s %9s | %-36s\n', 'Loads (FAST output)', 'BlPtchRt', 'BlMOoP', 'BlDEL', 'TwrDEL', ...
    'ratio PI/qLMPC: blPtchRt MOoP blDEL twrDEL');
fprintf('%-24s %9s %9s %9s %9s |\n', '', '(deg/s)', '(kNm)', '(kNm)', '(kNm)');
for k = 1:height(Tres)
    fprintf('%-24s %9.2f %9.0f %9.0f %9.0f | %5.2f %5.2f %5.2f %5.2f\n', Tres.Row{k}, ...
        Tres.BlPitchRateStd(k), Tres.BlMOoPStd(k), Tres.BlDEL(k), Tres.TwrDEL(k), ...
        Tres.ratio_BlPitchRateStd(k), Tres.ratio_BlMOoPStd(k), Tres.ratio_BlDEL(k), Tres.ratio_TwrDEL(k));
end
end

function writeTable(Tres, figDir, windName, windLabel, maxRate)
% writeTable writes the comparison table as csv and LaTeX table
writetable(Tres, fullfile(figDir, ['resultsFASTthesis_', windName, '.csv']), 'WriteRowNames', true);

fid = fopen(fullfile(figDir, ['resultsFASTthesis_', windName, '.tex']), 'w');
fprintf(fid, '\\begin{table}[h]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, ['\\caption[Controller comparison on FAST, %s.]{Controller comparison on FAST v8.16, %s, ', ...
    'pitch rate limit %g\\,\\si{\\degree/s}. Standard deviations (top) and blade pitch rate, ', ...
    'blade root out-of-plane moment standard deviation and damage equivalent loads of the blade root ', ...
    '(Woehler exponent 10, maximum of the blades) and tower base (exponent 4) (bottom); in parentheses ', ...
    'the ratio PI/qLMPC, values above 1 indicate improvement.}\n'], windLabel, windLabel, maxRate);
fprintf(fid, '\\label{tab:CtrlComparisonFAST%s}\n', windName);
% Control performance (std) and loads (blade pitch rate std, blade root
% out-of-plane moment std, DELs)
headers = {['{Controller} & {$\\omega_r$ (rpm)} & {$|P_{g,\\text{ref}}-P_g|$ (kW)} & ', ...
    '{Twr$_\\text{FA}$ (m/s$^2$)} & {$\\dot T_g$ (kNm/s)} & {$\\dot\\beta$ (\\si{\\degree/s})} \\\\\n'], ...
    ['{Controller} & {$\\dot\\beta_1$ (\\si{\\degree/s})} & {$M_{yb}$ (kNm)} & ', ...
    '{DEL$_b$ (kNm)} & {DEL$_t$ (kNm)} \\\\\n']};
metricSets = {{'RotSpdStd', 'PwrErrStd', 'TwrFAStd', 'GenTqRateStd', 'PitchRateStd'}, ...
    {'BlPitchRateStd', 'BlMOoPStd', 'BlDEL', 'TwrDEL'}};
fmtSets = {{'%.3f', '%.1f', '%.3f', '%.2f', '%.2f'}, {'%.2f', '%.0f', '%.0f', '%.0f'}};
for idxS = 1:2
    metricNames = metricSets{idxS};
    fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('c', 1, numel(metricNames)));
    fprintf(fid, '\\toprule\n');
    fprintf(fid, headers{idxS});
    fprintf(fid, '\\midrule\n');
    for k = 1:height(Tres)
        cells = cell(1, numel(metricNames));
        for m = 1:numel(metricNames)
            val = sprintf(fmtSets{idxS}{m}, Tres.(metricNames{m})(k));
            ratio = Tres.(['ratio_', metricNames{m}])(k);
            if k > 1 && isfinite(ratio)
                val = sprintf('%s (%.2f)', val, ratio);
            end
            cells{m} = val;
        end
        fprintf(fid, '%s & %s \\\\\n', Tres.Row{k}, strjoin(cells, ' & '));
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    if idxS == 1
        fprintf(fid, '\n\\vspace{2mm}\n\n');
    end
end
fprintf(fid, '\\end{table}\n');
fclose(fid);
end

function plotRuns(runs, ctrlNames, GenPwrRef, timeEval, figNo, figDir, windName, windLabel, fstFile, maxRate)
% plotRuns plots wind, generator torque, blade 1 pitch, tower fore-aft
% acceleration, generator power and blade 1 root out-of-plane moment of all
% runs in the evaluation window
DT = 0.008;
yAxCell = {'Wind (m/s)', 'GenTq (kNm)', 'Pitch_1 (°)', 'TwrAcc_{FA} (m/s^2)', 'GenPwr (MW)', 'M_{yb,1} (MNm)'};
[~, fstName] = fileparts(fstFile);
figure(figNo); clf;
set(gcf, 'Name', sprintf('PI vs. qLMPC designs (%s) | WECS: FAST v8.16 (%s.fst) | wind: %s | rate limit %g deg/s', ...
    strjoin(ctrlNames(2:end), ', '), fstName, windLabel, maxRate), 'NumberTitle', 'off');
tiledlayout(numel(yAxCell), 1, 'TileSpacing', 'Compact');
ax = gobjects(numel(yAxCell), 1);
for k = 1:numel(yAxCell)
    ax(k) = nexttile; hold on; grid on;
    ylabel(yAxCell{k});
end
lineStyles = {'-', '-.', '--', ':'};
for idxC = 1:numel(runs)
    [T, t] = runTable(runs{idxC}, DT);
    if isempty(T), continue, end
    idx = t >= timeEval(1) & t <= timeEval(2);
    ls = lineStyles{mod(idxC - 1, numel(lineStyles)) + 1};
    if idxC == 1
        plot(ax(1), t(idx), T.Wind(idx), 'k');
    end
    plot(ax(2), t(idx), T.GenTq(idx)/1000, ls, 'DisplayName', ctrlNames{idxC});
    plot(ax(4), t(idx), T.NcIMUTAxs(idx), ls);
    plot(ax(5), t(idx), T.GenPwr(idx)/1000, ls);
    % Blade 1 pitch and root moment from the FAST output (individual pitch)
    F = runs{idxC}.FASTOut;
    if ~isempty(F)
        idxF = F.Time >= timeEval(1) & F.Time <= timeEval(2);
        plot(ax(3), F.Time(idxF), F.BlPitch1(idxF), ls);
        plot(ax(6), F.Time(idxF), F.RootMOoP1(idxF)/1000, ls);
    else
        plot(ax(3), t(idx), T.BlPitch1(idx), ls);
        plot(ax(6), NaN, NaN, ls);
    end
end
tRef = (timeEval(1):DT:timeEval(2))';
tRef = tRef(tRef <= GenPwrRef(end, 1));
plot(ax(5), tRef, interp1(GenPwrRef(:, 1), GenPwrRef(:, 2), tRef, 'previous', 'extrap')/1000, 'k--');
title(ax(1), sprintf('FAST v8.16, %s, pitch rate limit %g°/s', windLabel, maxRate));
lgd = legend(ax(2));
lgd.Layout.Tile = 'east'; % outside the tiles, all axes keep the same width
xlabel(ax(end), 'Time (s)');
linkaxes(ax, 'x'); axis(ax, 'tight');
posDefault = get(0, 'DefaultFigurePosition');
set(gcf, 'position', [posDefault(1), posDefault(2) - posDefault(4)*0.7, posDefault(3)*1.4, posDefault(4)*1.7]);
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);
print(gcf, fullfile(figDir, ['cmpCtrlFASTthesis_', windName]), '-dpng');
print(gcf, fullfile(figDir, ['cmpCtrlFASTthesis_', windName]), '-depsc');
end
