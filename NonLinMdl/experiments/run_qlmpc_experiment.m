function M = run_qlmpc_experiment(name, plant, wind, params)
% run_qlmpc_experiment runs one closed-loop qLPV-MPC simulation with
% modified controller parameters and stores the evaluation.
%
% Inputs:
% - name:   experiment name (result file experiments/results/<name>.mat)
% - plant:  'FAST' (FAST v8.16, Mdl_FAST), 'SimFilt' (Mdl_BianchiOL_WindFiltered)
%           or 'SimPlain' (Mdl_BianchiOL)
% - wind:   'NTW18' (turbulent, 18 m/s mean) or 'Sweep' (4 to 25 m/s)
% - params: struct with base workspace variables that override the values
%           from initWorkspace, e.g. struct('qMPC', [...], 'fcOmegaCtrl', 0.25)
%
% The metrics are evaluated for 15 to 400 s (Simulink plants) or 60 to
% 600 s (FAST) and appended to
% experiments/results_log.txt.

if nargin < 4 || isempty(params)
    params = struct();
end

expDir = fileparts(mfilename('fullpath'));
workDir = fileparts(expDir);
mainDir = fileparts(workDir);
resDir = fullfile(expDir, 'results');
if ~isfolder(resDir), mkdir(resDir); end

cd(workDir);
evalin('base', 'run(''initWorkspace.m'')');
cache = fullfile(tempdir, 'qlmpc_experiments', name);
if ~isfolder(cache), mkdir(cache); end
Simulink.fileGenControl('set', 'CacheFolder', cache, 'CodeGenFolder', cache, 'createDir', true);

% Parameter overrides; ctrlVariant 'tower' selects the controller with the
% tower fore-aft model (MdlCtrl_qLPVMPCtower, qMPC with 14 elements),
% 'flap' the one with tower fore-aft and blade flap model
% (MdlCtrl_qLPVMPCflap, qMPC with 16 elements)
ctrlVariant = 'beta';
fn = fieldnames(params);
for k = 1:numel(fn)
    if strcmp(fn{k}, 'ctrlVariant')
        ctrlVariant = params.ctrlVariant;
    else
        assignin('base', fn{k}, params.(fn{k}));
    end
end

% Wind data (stop time, initial values and power reference)
switch wind
    case 'NTW18'
        S = load(fullfile(mainDir, 'dataIn', 'OutDataWind18NTW.mat'), 'OutTable');
        fst = 'FAST_NTW18.fst';
    case 'NTW18b' % second TurbSim seed (FAST only): TurbKlasseA/18mpsNewSeed.bts
        S = load(fullfile(mainDir, 'dataIn', 'OutDataWind18NTW.mat'), 'OutTable');
        fst = 'FAST_NTW18b.fst';
    case 'Sweep'
        S = load(fullfile(mainDir, 'dataIn', 'OutDataSweep.mat'), 'OutTable');
        fst = 'FAST.fst';
    otherwise
        error('Unknown wind case %s.', wind);
end
OutTable = S.OutTable;
assignin('base', 'OutTable', OutTable);

% Plant
fastFiles = {};
switch plant
    case 'FAST'
        mdl = 'test_FAST_qLPVMPCbeta';
        inDir = fullfile(workDir, 'FASTinputfiles');
        fstRun = fullfile(inDir, ['FAST_exp_' name '.fst']);   % own root name for parallel runs
        copyfile(fullfile(inDir, fst), fstRun);
        fastFiles = {fstRun};
        assignin('base', 'FAST_InputFileName', fstRun);
        assignin('base', 'TMax', OutTable.Time(end));
        assignin('base', 'OutList', cellstr(readlines(fullfile(inDir, 'OutList_FAST.txt'), 'EmptyLineRule', 'skip')));
    case {'SimFilt', 'SimPlain'}
        mdl = 'test_SimulinkMdl2_qLPVMPCbeta';
        load_system(mdl);
        blk = find_system(mdl, 'BlockType', 'ModelReference', 'Name', 'WECS Model');
        if strcmp(plant, 'SimFilt')
            set_param(blk{1}, 'ModelName', 'Mdl_BianchiOL_WindFiltered');
            assignin('base', 'step_thresshold', 1);
            assignin('base', 'kVeFilt', 0.25);
            assignin('base', 'kVeFilt_tau', 0.75);
        else
            set_param(blk{1}, 'ModelName', 'Mdl_BianchiOL');
        end
    otherwise
        error('Unknown plant %s.', plant);
end

if any(strcmp(ctrlVariant, {'tower', 'flap'}))
    if strcmp(plant, 'FAST')
        mdl = ['test_FAST_qLPVMPC' ctrlVariant];
    else
        close_system(mdl, 0);
        mdl = ['test_SimulinkMdl2_qLPVMPC' ctrlVariant];
        load_system(mdl);
        blk = find_system(mdl, 'BlockType', 'ModelReference', 'Name', 'WECS Model');
        if strcmp(plant, 'SimFilt')
            set_param(blk{1}, 'ModelName', 'Mdl_BianchiOL_WindFiltered');
        else
            set_param(blk{1}, 'ModelName', 'Mdl_BianchiOL');
        end
    end
    if ~evalin('base', 'exist(''fcTowerCtrl'', ''var'')')
        assignin('base', 'fcTowerCtrl', 0.02);
    end
    if ~evalin('base', 'exist(''towerSignCtrl'', ''var'')')
        assignin('base', 'towerSignCtrl', 1);
    end
end

% Simulation
M = struct('name', name, 'plant', plant, 'wind', wind, 'status', "ok", 'msg', "", 'ctrlVariant', ctrlVariant);
t0 = tic;
T = table();
errSig = [];
tictoc = [];
try
    load_system(mdl);
    % CaptureErrors: keep the logged data up to an abort (e.g. FAST divergence)
    out = sim(mdl, 'SrcWorkspace', 'base', 'CaptureErrors', 'on');
    if ~isempty(out.ErrorMessage)
        M.status = "failed";
        M.msg = string(out.ErrorMessage);
    end
    varnames = {'Wind', 'RotSpeed', 'GenPwr', 'GenTq', 'BlPitch1', 'NcIMUTAxs', 'NcIMUTAys'};
    D = out.get('OutDataTest');
    if ~isempty(D)
        T = array2table(D, 'VariableNames', varnames(1:size(D, 2)));
    end
    if ~isempty(out.get('err')), errSig = squeeze(out.err.Data); end
    if ~isempty(out.get('tictoc_LPVMPC')), tictoc = out.tictoc_LPVMPC.Data; end
catch ME
    M.status = "failed";
    M.msg = string(ME.message);
end
M.wallTime_s = toc(t0);
close_system(mdl, 0);

% Evaluation window: 15 to 400 s for the Simulink plant models (fast
% decaying transients), 60 to 600 s for FAST
DT = 0.008;
if strcmp(plant, 'FAST')
    M.tEval = [60 600];
else
    M.tEval = [15 400];
end
M.simulated_s = height(T) * DT;
if ~isempty(T)
    t = (0:height(T)-1)' * DT;
    idx = t >= M.tEval(1) & t <= M.tEval(2);
    dPitch = diff(T.BlPitch1(idx)) / DT;
    M.RotSpeed_mean = mean(T.RotSpeed(idx));
    M.RotSpeed_std  = std(T.RotSpeed(idx));
    M.RotSpeed_max  = max(T.RotSpeed(idx));
    M.GenPwr_mean   = mean(T.GenPwr(idx));
    M.GenPwr_std    = std(T.GenPwr(idx));
    M.TwrFA_std     = std(T.NcIMUTAxs(idx));
    M.PitchRate_std = std(dPitch);
    M.PitchRate_meanAbs = mean(abs(dPitch));
    M.GenTqRate_std = std(diff(T.GenTq(idx) / 1000) / DT);
    M.n_err         = nnz(errSig);
    M.tictoc_mean   = mean(tictoc);
    if M.status == "ok" && (M.simulated_s < OutTable.Time(end) - 1 || M.RotSpeed_max > 15)
        M.status = "unstable";
    end
end
T5 = T(1:5:end, :);   % 0.04 s resolution
save(fullfile(resDir, [name '.mat']), 'M', 'params', 'T5');

% Log line
logFile = fullfile(expDir, 'results_log.txt');
if ~isfile(logFile)
    fid = fopen(logFile, 'w');
    fprintf(fid, ['name\tplant\twind\tstatus\tsimulated_s\tRotSpeed_mean\tRotSpeed_std\tRotSpeed_max\tGenPwr_mean\tGenPwr_std\t' ...
        'TwrFA_std\tPitchRate_std\tGenTqRate_std\tn_err\tparams\n']);
    fclose(fid);
end
v = @(f) getfieldNaN(M, f);
fid = fopen(logFile, 'a');
fprintf(fid, '%s\t%s\t%s\t%s\t%.1f\t%.3f\t%.3f\t%.3f\t%.0f\t%.1f\t%.4f\t%.3f\t%.3f\t%d\t%s\n', name, plant, wind, M.status, ...
    M.simulated_s, v('RotSpeed_mean'), v('RotSpeed_std'), v('RotSpeed_max'), v('GenPwr_mean'), v('GenPwr_std'), ...
    v('TwrFA_std'), v('PitchRate_std'), v('GenTqRate_std'), v('n_err'), jsonencode(params));
fclose(fid);
disp(M);

% Remove the FAST copies and outputs of this run
for k = 1:numel(fastFiles)
    [d, r] = fileparts(fastFiles{k});
    delete(fullfile(d, [r '.*']));
end
end

function x = getfieldNaN(S, f)
if isfield(S, f), x = S.(f); else, x = NaN; end
end
