%% Script to compare the PI and the most promising qLPV-MPC designs on FAST
% Closed-loop runs of the baseline PI controller and the stored qLPV-MPC
% designs (subfunctionsMPC/qlmpcDesign.m) with FAST v8.16 as the WECS model
% for the turbulent 18 m/s wind (FAST_NTW18.fst, TurbSim) and the wind sweep
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
% Metrics (evaluation window 60 to 600 s, limited by the data length): std
% of the rotor speed, the power tracking error |P_g,ref - P_g|, the tower
% fore-aft acceleration, the generator torque rate and the pitch rate; the
% ratio sigma_PI/sigma_qLMPC is above 1 if the qLMPC design is better.
clc; clear; close all;

Simulink.data.dictionary.closeAll('-discard');

%% Set path to initialization script initWorkspace and run it
addpath(genpath('NonLinMdl'));
addpath(genpath('LinMdl'));

initWorkspace;

%% Settings
loadData = 1; % load the stored simulation results if available
maxRate = 13; % pitch rate limit (deg/s) for PI and qLMPC

% qLPV-MPC designs (qlmpcDesign.m, results in NonLinMdl/experiments/EXPERIMENTS.md)
% - 'omega': rotor speed tracking, lowest pitch activity and tower acceleration
% - 'flapNotch': tower/blade flap model with notch filters, best speed and power
% - 'flapNotchSmooth12': DTc = 0.12 s, smooth pitch reference, balanced
% - 'flapNotchSmooth12Tower': as above with the tower velocity weight
designNames = {'omega', 'flapNotch', 'flapNotchSmooth12', 'flapNotchSmooth12Tower'};

% Wind cases: name, FASTTool data (initial values, stop time), FAST input file
windCases = {'NTW18', 'OutDataWind18NTW.mat', 'FAST_NTW18.fst', '18 m/s turbulent';
    'Sweep', 'OutDataSweep.mat', 'FAST.fst', 'sweep'};

timeEval = [60, 600]; % evaluation window (s), the first 60 s are disregarded
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

    % Run or load the simulations
    runs = cell(size(ctrlNames));
    for idxC = 1:numel(ctrlNames)
        matFile = fullfile(dataDirOut, sprintf('OutTable_%s_%s.mat', windName, ctrlNames{idxC}));
        runs{idxC} = getFASTRun(matFile, loadData, ctrlNames{idxC}, OutTable, fstFile, maxRate);
    end

    % Power reference of the qLMPC (depends on the measured wind only, the
    % FAST wind field is the same for all controllers); used for all runs
    idxRef = find(cellfun(@(R) ~isempty(R.GenPwrRef), runs), 1);
    if isempty(idxRef)
        error('No power reference logged for wind case %s.', windName);
    end
    GenPwrRef = runs{idxRef}.GenPwrRef;

    % Metrics
    resultsFAST.(windName) = evaluateRuns(runs, ctrlNames, GenPwrRef, timeEval);
    printTable(resultsFAST.(windName), windCases{idxW, 4}, maxRate);
    writeTable(resultsFAST.(windName), figDir, windName, windCases{idxW, 4}, maxRate);

    % Time series
    plotRuns(runs, ctrlNames, GenPwrRef, timeEval, figNo, figDir, windName, windCases{idxW, 4}, fstFile, maxRate);
    figNo = figNo + 1;
end
save(fullfile(dataDirOut, 'resultsFASTthesis.mat'), 'resultsFAST', 'designNames', 'windCases', 'maxRate', 'timeEval');

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
% struct, empty for the PI), status, msg, maxRate, wallTime_s
if loadData && isfile(matFile)
    tmp = load(matFile, 'R');
    if ~isempty(tmp.R.OutTableRun) % runs without any data are repeated
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
    'maxRate', maxRate, 'fstFile', fstFile, 'OutTableRun', table(), 'GenPwrRef', [], 'tictoc', []);
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
delete(fullfile(fileparts(fstRun), 'FAST_exp_thesis.*'));
if R.status ~= "ok"
    warning('%s: %s', ctrlName, R.msg);
end
save(matFile, 'R');
fprintf('Saved %s (%.0f s)\n', matFile, R.wallTime_s);
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
    tCpuMean, simulated] = deal(nan(nC, 1));
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
end
% Numerically constant signals (constant torque above rated) count as zero
GenTqRateStd(GenTqRateStd < 1e-6) = 0;
Tres = table(RotSpdStd, RotSpdMax, PwrErrStd, TwrFAStd, GenTqRateStd, PitchRateStd, ...
    tCpuMean, simulated, status, 'RowNames', ctrlNames(:));
% Std ratios PI/qLMPC (above 1: qLMPC better)
metricNames = {'RotSpdStd', 'PwrErrStd', 'TwrFAStd', 'GenTqRateStd', 'PitchRateStd'};
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
end

function writeTable(Tres, figDir, windName, windLabel, maxRate)
% writeTable writes the comparison table as csv and LaTeX table
writetable(Tres, fullfile(figDir, ['resultsFASTthesis_', windName, '.csv']), 'WriteRowNames', true);

fid = fopen(fullfile(figDir, ['resultsFASTthesis_', windName, '.tex']), 'w');
fprintf(fid, '\\begin{table}[h]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, ['\\caption[Controller comparison on FAST, %s.]{Controller comparison on FAST v8.16, %s, ', ...
    'pitch rate limit %g\\,\\si{\\degree/s}. Standard deviations; in parentheses the ratio ', ...
    '$\\sigma_\\text{PI}/\\sigma_\\text{qLMPC}$, values above 1 indicate improvement.}\n'], windLabel, windLabel, maxRate);
fprintf(fid, '\\label{tab:CtrlComparisonFAST%s}\n', windName);
fprintf(fid, '\\begin{tabular}{lccccc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, ['{Controller} & {$\\omega_r$ (rpm)} & {$|P_{g,\\text{ref}}-P_g|$ (kW)} & ', ...
    '{Twr$_\\text{FA}$ (m/s$^2$)} & {$\\dot T_g$ (kNm/s)} & {$\\dot\\beta$ (\\si{\\degree/s})} \\\\\n']);
fprintf(fid, '\\midrule\n');
metricNames = {'RotSpdStd', 'PwrErrStd', 'TwrFAStd', 'GenTqRateStd', 'PitchRateStd'};
fmt = {'%.3f', '%.1f', '%.3f', '%.2f', '%.2f'};
for k = 1:height(Tres)
    cells = cell(1, numel(metricNames));
    for m = 1:numel(metricNames)
        val = sprintf(fmt{m}, Tres.(metricNames{m})(k));
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
fprintf(fid, '\\end{table}\n');
fclose(fid);
end

function plotRuns(runs, ctrlNames, GenPwrRef, timeEval, figNo, figDir, windName, windLabel, fstFile, maxRate)
% plotRuns plots wind, generator torque, pitch, tower fore-aft acceleration
% and generator power of all runs in the evaluation window
DT = 0.008;
yAxCell = {'Wind (m/s)', 'GenTq (kNm)', 'Pitch (°)', 'TwrAcc_{FA} (m/s^2)', 'GenPwr (MW)'};
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
lineStyles = {'-', '-.', '--', ':', '-.'};
for idxC = 1:numel(runs)
    [T, t] = runTable(runs{idxC}, DT);
    if isempty(T), continue, end
    idx = t >= timeEval(1) & t <= timeEval(2);
    ls = lineStyles{min(idxC, numel(lineStyles))};
    if idxC == 1
        plot(ax(1), t(idx), T.Wind(idx), 'k');
    end
    plot(ax(2), t(idx), T.GenTq(idx)/1000, ls, 'DisplayName', ctrlNames{idxC});
    plot(ax(3), t(idx), T.BlPitch1(idx), ls);
    plot(ax(4), t(idx), T.NcIMUTAxs(idx), ls);
    plot(ax(5), t(idx), T.GenPwr(idx)/1000, ls);
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
