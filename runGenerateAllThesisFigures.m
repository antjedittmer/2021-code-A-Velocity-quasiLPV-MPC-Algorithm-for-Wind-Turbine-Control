%% Script to generate all figures of the thesis created with MATLAB
% Calls the figure scripts of the three repositories next to this one
% (same parent folder):
%   2021-code-A-Velocity-quasiLPV-MPC-Algorithm-for-Wind-Turbine-Control
%     - runGeneratePicsPaper: Bode plots, model and controller comparison
%       on Mdl1/Mdl2/FAST, pitch rate constraints (figDir, figDirConstr1)
%     - runGeneratePicsFASTthesis: PI and qLMPC on FAST (figDirFASTthesis)
%     - plotParetoVerif2: Pareto and bar plot qLMPC/IPC/AALC (figDirPareto)
%   Combined-IPC-TEFC-for-Load-Reduction-of-Wind-Turbines (Main/plotFunktionen)
%     - mainPlotPaper: Reynolds number, lift and drag coefficients, TEF
%       open loop (fig1, fig2)
%     - plotIPC_RateLimitComparison: IPC time series, IPC rate limit PSD
%       and AALC PSDs (figIPCThesis)
%   2023-code-IFAC-Koopman-Model-Predictive-Control-for-Wind-Farm
%     - mainGeneratePlotThesis: Koopman WRC/AIC (figures_IFAC23)
%     - mainGeneratePlotThesisAdaptiveAIC: adaptive AIC (figures_CDC22)
% The scripts load stored simulation results if available (see the
% scripts), so a call mainly replots.
%
% Afterwards all png and eps files of the figure folders are copied into
% the folder AllFiguresThesis next to this repository, one subfolder per
% thesis figure folder (thesis_DLR_20260630/figures/<subfolder>), e.g.
% AllFiguresThesis/Figures_ADi. Existing files there are overwritten.
%
% Each script runs in the base workspace (Simulink models read their
% variables from there) and may clear it; the state of this script is kept
% in the local function main. A failing script is reported and skipped;
% the MATLAB path and the current folder are restored after each script.
%
% Thesis figures without a generating script in these repositories (copy
% them from the thesis folder if needed): figDirECC21 (Cp, Cq, Ct; plotOn
% in initModel5MWNREL), figDirETC, figures_IFAC23/01FigPwrVsYaw_tex,
% 03FigAngle0_20_tex, gridmesh2_tex, figuresBlockDiagram_IFAC23 and the
% pdf/jpg drawings.

main();

function main()
%% Steps: {repository, folder in the repository, script or function call, description}
% Set the last column to false to skip a step (e.g. the long Koopman runs).
mainDir = fileparts(mfilename('fullpath'));
parentDir = fileparts(mainDir);
velDir = mainDir;
ipcDir = fullfile(parentDir, 'Combined-IPC-TEFC-for-Load-Reduction-of-Wind-Turbines', 'Main', 'plotFunktionen');
koopDir = fullfile(parentDir, '2023-code-IFAC-Koopman-Model-Predictive-Control-for-Wind-Farm');

steps = {
    velDir,  'runGeneratePicsPaper.m',                 'Mdl1/Mdl2/FAST, controller comparison, pitch rate', true
    velDir,  'runGeneratePicsFASTthesis.m',            'PI and qLMPC on FAST', true
    velDir,  'plotParetoVerif2()',                     'Pareto and bar plot', true
    ipcDir,  'mainPlotPaper.m',                        'Re, airfoil coefficients, TEF open loop', true
    ipcDir,  'plotIPC_RateLimitComparison()',          'IPC time series, IPC and AALC PSDs', true
    koopDir, 'mainGeneratePlotThesis.m',               'Koopman WRC/AIC', true
    koopDir, 'mainGeneratePlotThesisAdaptiveAIC.m',    'adaptive AIC', true};

%% Figure folders: {source folder, subfolder in AllFiguresThesis (thesis folder name)}
figFolders = {
    fullfile(velDir, 'figDir'),              'Figures_ADi'
    fullfile(velDir, 'figDirConstr1'),       'figDirConstr1'
    fullfile(velDir, 'figDirFASTthesis'),    'figDirFASTthesis'
    fullfile(velDir, 'figDirPareto'),        'figures_WESC25'
    fullfile(ipcDir, 'fig1'),                'Images_Torque24'
    fullfile(ipcDir, 'fig2'),                'Images_Torque24'
    fullfile(ipcDir, 'figIPCThesis'),        'figIPCThesis'
    fullfile(koopDir, 'figures_IFAC23'),     'figures_IFAC23'
    fullfile(koopDir, 'figures_CDC22'),      'figures_CDC22'};
outDir = fullfile(parentDir, 'AllFiguresThesis');

%% Run the figure scripts
status = strings(size(steps, 1), 1);
for idx = 1:size(steps, 1)
    [repoDir, cmd, descr, doRun] = steps{idx, :};
    if ~doRun
        status(idx) = "skipped";
        continue
    end
    fprintf('\n===== %d/%d %s (%s) =====\n', idx, size(steps, 1), cmd, descr);
    status(idx) = runStep(repoDir, cmd);
end

%% Copy the png and eps files
nCopied = 0;
for idx = 1:size(figFolders, 1)
    [srcDir, subDir] = figFolders{idx, :};
    if ~isfolder(srcDir)
        warning('Figure folder %s not found', srcDir);
        continue
    end
    dstDir = fullfile(outDir, subDir);
    if ~isfolder(dstDir)
        mkdir(dstDir);
    end
    files = [dir(fullfile(srcDir, '*.png')); dir(fullfile(srcDir, '*.eps'))];
    for idxF = 1:numel(files)
        copyfile(fullfile(files(idxF).folder, files(idxF).name), fullfile(dstDir, files(idxF).name));
    end
    nCopied = nCopied + numel(files);
end

%% Summary
fprintf('\n===== Summary =====\n');
for idx = 1:size(steps, 1)
    fprintf('%-40s %s\n', steps{idx, 2}, status(idx));
end
fprintf('%d png and eps files copied to %s\n', nCopied, outDir);
end

function status = runStep(repoDir, cmd)
% Runs cmd (script file or function call) in repoDir in the base workspace;
% restores the MATLAB path and the current folder afterwards
oldPath = path;
oldDir = pwd;
cleanup = onCleanup(@() restoreState(oldPath, oldDir));
try
    cd(repoDir);
    if endsWith(cmd, '.m')
        evalin('base', sprintf('run(''%s'');', fullfile(repoDir, cmd)));
    else
        evalin('base', [cmd, ';']);
    end
    status = "ok";
catch ME
    fprintf(2, 'Error in %s: %s\n', cmd, ME.message);
    status = "failed: " + string(ME.message);
end
close all
end

function restoreState(oldPath, oldDir)
path(oldPath);
cd(oldDir);
end
