%% Script to generate all figures of the thesis created with MATLAB

% Titles of the figures (with diagnostics such as VAF values, gains or
% colour keys) per step of the step list in main, in this order:
%  1 runGeneratePicsPaper           2 runGeneratePicsFASTthesis
%  3 plotParetoVerif2               4 c_T, c_P, c_Q maps (no titles)
%  5 mainPlotPaper (IPC/TEFC)       6 plotIPC_RateLimitComparison
%  7 runWFSimAnimationPlots (none)  8 mainGeneratePlotThesis
%  9 mainGeneratePlotThesisAdaptiveAIC
% 10 WFSim grids and program structure (no titles)
% 11 plotScholarStatisticsWT       12 plotScholarStatisticsWF
% false: no titles (the title calls of the plot functions are skipped),
% true: titles as before. Scripts run on their own show the titles.
showTitles = [false, false, false, false, false, false, false, false, false, false, true, true];

% Font size of the time series figures of all repositories (model and
% controller comparisons, pitch rate, FAST, IPC, Koopman open and closed
% loop; getTimeSeriesFontSize). Empty: default of the scripts (13); a
% number overrides it for this call. Scripts run on their own use 13.
timeSeriesFontSize = [];

% Font size of the time axis labels ('Time (s)') of the time series figures
% (setTimeLabelFontSize), larger than the other texts. Empty: default of
% the scripts (14); a number overrides it for this call.
timeLabelFontSize = [];

% Calls the figure scripts of the repositories next to this one (same
% parent folder):
%   2021-code-A-Velocity-quasiLPV-MPC-Algorithm-for-Wind-Turbine-Control
%     - runGeneratePicsPaper: Bode plots, model and controller comparison
%       on Mdl1/Mdl2/FAST, pitch rate constraints (figDir, figDirConstr1)
%     - runGeneratePicsFASTthesis: PI and qLMPC on FAST (figDirFASTthesis)
%     - plotParetoVerif2: Pareto and bar plot qLMPC/IPC/AALC (figDirPareto)
%     - initModel5MWNREL: c_T, c_P and c_Q maps of the original NREL 5 MW
%       look-up tables (model version 'paper', no Cq correction;
%       figDirECC21)
%   Combined-IPC-TEFC-for-Load-Reduction-of-Wind-Turbines (Main/plotFunktionen)
%     - mainPlotPaper: Reynolds number, lift and drag coefficients, TEF
%       open loop (fig1, fig2)
%     - plotIPC_RateLimitComparison: IPC time series, IPC rate limit PSD
%       and AALC PSDs (figIPCThesis)
%   2023-code-IFAC-Koopman-Model-Predictive-Control-for-Wind-Farm
%     - runWFSimAnimationPlots (as in mainGeneratePlots): power over yaw
%       and wake centre lines (figPaper; PowerVsYaw and
%       gamma0_20_AngleCenter are the thesis figures 01FigPwrVsYaw_tex and
%       03FigAngle0_20_tex)
%     - mainGeneratePlotThesis: Koopman WRC/AIC (figures_IFAC23)
%     - mainGeneratePlotThesisAdaptiveAIC: adaptive AIC (figures_CDC22)
%     - WFSimCode/bin/analysis: plotMeshDiss_zoombox, grids of the two and
%       nine turbine farms (dirFigWFSimLayout; thesis figures_IFAC23 and
%       figDirETC), plotFlowchartInitCtrlPlant, WFSim program structure
%       (figuresBlockDiagram_IFAC23)
%   Scholar (https://github.com/antjedittmer/Scholar, MATLABscripts4plots)
%     - plotScholarStatisticsWT, plotScholarStatisticsWF: Google Scholar
%       entries per year (dirFigures)
% The scripts load stored simulation results if available (see the
% scripts), so a call mainly replots.
%
% Afterwards all png and eps files of the figure folders are copied into
% two folders next to this repository, depending on whether the thesis
% (thesisDir/dissertation_ADi.tex and the files it includes or inputs)
% links the figure (png, eps or pdf of the same name):
%   AllFiguresThesisUsed/<subfolder>, AllFiguresThesisUnused/<subfolder>
% <subfolder> is the thesis figure folder (thesis figures/<subfolder>),
% e.g. Figures_ADi. Both folders are deleted and created again on each
% call. Thesis figures (eps, png) that none of the scripts creates are
% listed at the end.
%
% Titles: the vector showTitles in the first lines switches the titles of
% the figures per step (default: off, on for the Google Scholar plots).
%
% Each script runs in the base workspace (Simulink models read their
% variables from there) and may clear it; the state of this script is kept
% in the local function main. A failing script is reported and skipped;
% the MATLAB path and the current folder are restored after each script.

main(showTitles, timeSeriesFontSize, timeLabelFontSize);

function main(showTitles, timeSeriesFontSize, timeLabelFontSize)
%% Folders
mainDir = fileparts(mfilename('fullpath'));
parentDir = fileparts(mainDir);
velDir = mainDir;
ipcDir = fullfile(parentDir, 'Combined-IPC-TEFC-for-Load-Reduction-of-Wind-Turbines', 'Main', 'plotFunktionen');
if exist(ipcDir,'dir') ~= 7 % ipcDir is on a different location on my PC
    studentDir = fullfile(parentDir,'Studenten','2023_BernhardOsterwinter','code');
    ipcDir = strrep(ipcDir,parentDir,studentDir);
    if exist(ipcDir,'dir') ~= 7,error('ipcdir does not exist. Please check the location'); end
end   
koopDir = fullfile(parentDir, '2023-code-IFAC-Koopman-Model-Predictive-Control-for-Wind-Farm');
if exist(koopDir,'dir') ~= 7, error('koopDir does not exist. Please check the location'); end
scholarDir = fullfile(parentDir, 'Scholar', 'MATLABscripts4plots'); % MATLABprojects\Scholar
if exist(scholarDir,'dir') ~= 7 % scholarDiris on a different location on my PC
    scholarDirParent = fullfile(parentDir,'MATLABprojects');
    scholarDir = strrep(scholarDir,parentDir,scholarDirParent);
    if exist(scholarDir,'dir') ~= 7,error('scholarDir does not exist. Please check the location'); end
end 
thesisDir = fullfile(parentDir, 'thesis_DLR_20260630');
if exist(thesisDir,'dir') ~= 7 % thesisDir is on a different location on my PC
    thesisDirParent = fullfile(fileparts(parentDir),'Documents','Dissertation');
    thesisDir = strrep(thesisDir,parentDir, thesisDirParent);
    if exist(thesisDir,'dir') ~= 7,error('scholarDir does not exist. Please check the location'); end
end 

usedDir = fullfile(parentDir, 'AllFiguresThesisUsed'); % Outputdir one
unusedDir = fullfile(parentDir, 'AllFiguresThesisUnused'); % Outputdir one

%% Steps: {folder, script, function call or function handle, description, run}
% Set the last column to false to skip a step (e.g. the long Koopman runs).
steps = {
    velDir,     'runGeneratePicsPaper.m',               'Mdl1/Mdl2/FAST, controller comparison, pitch rate', true
    velDir,     'runGeneratePicsFASTthesis.m',          'PI and qLMPC on FAST', true
    velDir,     'plotParetoVerif2()',                   'Pareto and bar plot', true
    velDir,     @() plotCoefficientMaps(velDir),        'c_T, c_P, c_Q maps (initModel5MWNREL)', true
    ipcDir,     'mainPlotPaper.m',                      'Re, airfoil coefficients, TEF open loop', true
    ipcDir,     'plotIPC_RateLimitComparison()',        'IPC time series, IPC and AALC PSDs', true
    koopDir,    fullfile('WFSimCode', 'bin', 'analysis', 'runWFSimAnimationPlots.m'), 'power over yaw, wake centre lines', true
    koopDir,    'mainGeneratePlotThesis.m',             'Koopman WRC/AIC', true
    koopDir,    'mainGeneratePlotThesisAdaptiveAIC.m',  'adaptive AIC', true
    koopDir,    @() plotWFSimStructure(koopDir),        'WFSim grids and program structure', true
    scholarDir, 'plotScholarStatisticsWT.m',            'Google Scholar entries wind turbine control', true
    scholarDir, 'plotScholarStatisticsWF.m',            'Google Scholar entries wind farm control', true};

%% Figure folders: {source folder, thesis subfolder}
figFolders = {
    fullfile(velDir, 'figDir'),                     'Figures_ADi'
    fullfile(velDir, 'figDirConstr1'),              'figDirConstr1'
    fullfile(velDir, 'figDirFASTthesis'),           'figDirFASTthesis'
    fullfile(velDir, 'figDirPareto'),               'figures_WESC25'
    fullfile(velDir, 'figDirECC21'),                'figDirECC21'
    fullfile(ipcDir, 'fig1'),                       'Images_Torque24'
    fullfile(ipcDir, 'fig2'),                       'Images_Torque24'
    fullfile(ipcDir, 'figIPCThesis'),               'figIPCThesis'
    fullfile(koopDir, 'figPaper'),                  'figPaper'
    fullfile(koopDir, 'figures_IFAC23'),            'figures_IFAC23'
    fullfile(koopDir, 'figures_CDC22'),             'figures_CDC22'
    fullfile(koopDir, 'WFSimCode', 'bin', 'analysis', 'figuresBlockDiagram_IFAC23'), 'figuresBlockDiagram_IFAC23'
    fullfile(fileparts(scholarDir), 'dirFigures'),  'Figures_ADi'};
% Figures with another name in the thesis: {source without extension,
% thesis subfolder, thesis name}
figRenamed = {
    fullfile(koopDir, 'figPaper', 'PowerVsYaw'),            'figures_IFAC23', '01FigPwrVsYaw_tex'
    fullfile(koopDir, 'figPaper', 'gamma0_20_AngleCenter'), 'figures_IFAC23', '03FigAngle0_20_tex'
    fullfile(koopDir, 'WFSimCode', 'bin', 'analysis', 'dirFigWFSimLayout', 'gridmesh2_tex'), 'figures_IFAC23', 'gridmesh2_tex'
    fullfile(koopDir, 'WFSimCode', 'bin', 'analysis', 'dirFigWFSimLayout', 'gridmesh9_tex'), 'figDirETC', 'gridmesh9_tex'};

%% Run the figure scripts
nSteps = size(steps, 1);
if numel(showTitles) ~= nSteps
    error('showTitles must have %d elements, one per step.', nSteps);
end
status = strings(nSteps, 1);
for idx = 1:nSteps
    [repoDir, cmd, descr, doRun] = steps{idx, :};
    if ~doRun
        status(idx) = "skipped";
        continue
    end
    fprintf('\n===== %d/%d %s =====\n', idx, nSteps, descr);
    status(idx) = runStep(repoDir, cmd, showTitles(idx), timeSeriesFontSize, timeLabelFontSize);
end

%% Figures linked in the thesis (subfolder/name, lower case, no extension)
linked = thesisFigures(fullfile(thesisDir, 'dissertation_ADi.tex'));

%% Copy the png and eps files into the used and unused folder
for outDir = {usedDir, unusedDir}
    if isfolder(outDir{1})
        rmdir(outDir{1}, 's');
    end
    mkdir(outDir{1});
end
src = {}; dst = {}; % source files, target subfolder/name
for idx = 1:size(figFolders, 1)
    [srcDir, subDir] = figFolders{idx, :};
    if ~isfolder(srcDir)
        warning('Figure folder %s not found', srcDir);
        continue
    end
    files = [dir(fullfile(srcDir, '*.png')); dir(fullfile(srcDir, '*.eps'))];
    for idxF = 1:numel(files)
        src{end+1} = fullfile(files(idxF).folder, files(idxF).name); %#ok<AGROW>
        dst{end+1} = [subDir, '/', files(idxF).name]; %#ok<AGROW>
    end
end
for idx = 1:size(figRenamed, 1)
    for ext = {'.png', '.eps'}
        if isfile([figRenamed{idx, 1}, ext{1}])
            src{end+1} = [figRenamed{idx, 1}, ext{1}]; %#ok<AGROW>
            dst{end+1} = [figRenamed{idx, 2}, '/', figRenamed{idx, 3}, ext{1}]; %#ok<AGROW>
        end
    end
end
isUsed = false(size(dst));
for idx = 1:numel(dst)
    [subDir, name] = fileparts(dst{idx});
    isUsed(idx) = ismember(lower([subDir, '/', name]), linked);
    if isUsed(idx), outDir = usedDir; else, outDir = unusedDir; end
    if ~isfolder(fullfile(outDir, subDir))
        mkdir(fullfile(outDir, subDir));
    end
    copyfile(src{idx}, fullfile(outDir, subDir, [name, dst{idx}(end-3:end)]));
end

%% Summary
fprintf('\n===== Summary =====\n');
for idx = 1:nSteps
    fprintf('%-50s %s\n', steps{idx, 3}, status(idx));
end
fprintf('%d png and eps files linked in the thesis copied to %s\n', nnz(isUsed), usedDir);
fprintf('%d png and eps files not linked in the thesis copied to %s\n', nnz(~isUsed), unusedDir);
created = unique(cellfun(@(f) lower(f(1:end-4)), dst, 'UniformOutput', false));
linkedEpsPng = thesisFigures(fullfile(thesisDir, 'dissertation_ADi.tex'), {'.eps', '.png'});
missing = setdiff(linkedEpsPng, created);
if ~isempty(missing)
    fprintf('Thesis figures (eps, png) not created by the scripts:\n');
    fprintf('  %s\n', missing{:});
end
end

function status = runStep(repoDir, cmd, showTitle, timeSeriesFontSize, timeLabelFontSize)
% Runs cmd in repoDir: a script file (path relative to repoDir) or a
% function call in the base workspace, or a function handle; restores the
% MATLAB path and the current folder afterwards. showTitle: titles of the
% figures (root application data ShowThesisTitle, read by the plot
% functions, removed afterwards); timeSeriesFontSize: font size of the time
% series figures (root application data TimeSeriesFontSize, read by
% getTimeSeriesFontSize, removed afterwards; empty: default 13);
% timeLabelFontSize: font size of the time axis labels (root application
% data TimeLabelFontSize, read by setTimeLabelFontSize, removed afterwards;
% empty: default 14)
oldPath = path;
oldDir = pwd;
setappdata(groot, 'ShowThesisTitle', logical(showTitle));
if ~isempty(timeSeriesFontSize)
    setappdata(groot, 'TimeSeriesFontSize', timeSeriesFontSize);
end
if ~isempty(timeLabelFontSize)
    setappdata(groot, 'TimeLabelFontSize', timeLabelFontSize);
end
cleanup = onCleanup(@() restoreState(oldPath, oldDir));
try
    if ~isfolder(repoDir)
        error('Folder %s not found', repoDir);
    end
    cd(repoDir);
    if isa(cmd, 'function_handle')
        cmd();
    elseif endsWith(cmd, '.m')
        evalin('base', sprintf('run(''%s'');', fullfile(repoDir, cmd)));
    else
        evalin('base', [cmd, ';']);
    end
    status = "ok";
catch ME
    fprintf(2, 'Error: %s\n', ME.message);
    status = "failed: " + string(ME.message);
end
close all
end

function restoreState(oldPath, oldDir)
path(oldPath);
cd(oldDir);
if isappdata(groot, 'ShowThesisTitle')
    rmappdata(groot, 'ShowThesisTitle');
end
if isappdata(groot, 'TimeSeriesFontSize')
    rmappdata(groot, 'TimeSeriesFontSize');
end
if isappdata(groot, 'TimeLabelFontSize')
    rmappdata(groot, 'TimeLabelFontSize');
end
end

function plotCoefficientMaps(velDir)
% c_T, c_P and c_Q maps of initModel5MWNREL in figDirECC21 with the model
% version 'paper' (original look-up tables, as in the thesis); the base
% workspace variable wecsModelVersion is restored afterwards
addpath(genpath(fullfile(velDir, 'LinMdl')));
addpath(genpath(fullfile(velDir, 'NonLinMdl')));
hasVersion = evalin('base', 'exist(''wecsModelVersion'', ''var'')');
oldVersion = '';
if hasVersion
    oldVersion = evalin('base', 'wecsModelVersion');
end
assignin('base', 'wecsModelVersion', 'paper');
cleanup = onCleanup(@() restoreVersion(hasVersion, oldVersion));
figDir = fullfile(velDir, 'figDirECC21');
if ~isfolder(figDir)
    mkdir(figDir);
end
initModel5MWNREL(1, '', '', '', '', figDir);
end

function plotWFSimStructure(koopDir)
% Grids of the two and nine turbine farms with zoomed turbine regions
% (plotMeshDiss_zoombox, dirFigWFSimLayout) and the WFSim program structure
% diagrams (plotFlowchartInitCtrlPlant, figuresBlockDiagram_IFAC23) in
% WFSimCode/bin/analysis; runStep restores the current folder
cd(fullfile(koopDir, 'WFSimCode', 'bin', 'analysis'));
plotMeshDiss_zoombox(2, false);
plotMeshDiss_zoombox(9, false);
plotFlowchartInitCtrlPlant(1);
plotFlowchartInitCtrlPlant(2);
end

function restoreVersion(hasVersion, oldVersion)
if hasVersion
    assignin('base', 'wecsModelVersion', oldVersion);
else
    evalin('base', 'clear wecsModelVersion');
end
end

function linked = thesisFigures(texFile, exts)
% Figures of \includegraphics in texFile and the files it includes or
% inputs (recursively, without comments) as lower case 'subfolder/name'
% (no extension, without the leading 'figures/'); exts: only these
% extensions (default: all)
if nargin < 2, exts = {}; end
linked = {};
if ~isfile(texFile)
    warning('Thesis file %s not found', texFile);
    return
end
rootDir = fileparts(texFile);
todo = {texFile};
done = {};
while ~isempty(todo)
    f = todo{1}; todo(1) = [];
    if ismember(f, done) || ~isfile(f), continue, end
    done{end+1} = f; %#ok<AGROW>
    txt = regexprep(fileread(f), '(?<!\\)%[^\n]*', ''); % comments
    inc = regexp(txt, '\\(?:include|input)\s*\{([^}]+)\}', 'tokens');
    for k = 1:numel(inc)
        g = fullfile(rootDir, strtrim(inc{k}{1}));
        if ~endsWith(g, '.tex'), g = [g, '.tex']; end %#ok<AGROW>
        todo{end+1} = g; %#ok<AGROW>
    end
    figs = regexp(txt, '\\includegraphics\s*(?:\[[^\]]*\])?\s*\{([^}]+)\}', 'tokens');
    for k = 1:numel(figs)
        [p, name, ext] = fileparts(strrep(strtrim(figs{k}{1}), '\', '/'));
        if ~isempty(exts) && ~ismember(lower(ext), exts), continue, end
        p = regexprep(p, '^figures/', '', 'ignorecase'); % not Figures_ADi
        linked{end+1} = lower([p, '/', name]); %#ok<AGROW>
    end
end
linked = unique(linked);
end
