%% Generate the FAST reference dataIn/OutDataWind18NTW.mat (18 m/s turbulent wind)
% FAST v8.16 (S-function of this repository) with the PI baseline
% (test_FAST_Baseline) in the 18 m/s normal turbulence wind field
% FASTinputfiles/wind18_NTM_B.bts: IEC 61400-1 Ed3, class B, NTM, power law
% shear 0.2, 660 s, 20x20 grid of 126 m (parameters of wind.inp, generated
% with TurbSim v3.2.1 from wind18_NTM_B.inp). FAST input: FAST_NTM18B.fst.
%
% The table contains the FAST output channels of OutList_FAST.txt in FAST
% units (as the earlier FASTTool file), DT = 0.008 s, 0 to 660 s.
% The wind of the test model (OutTable.Wind1VelX) only sets the initial
% values (operating point lookup tables), so the run is repeated with the
% hub height wind of the first run.

initWorkspace;
dirWork = fileparts(which('initWorkspace'));
inDir = fullfile(dirWork, 'FASTinputfiles');
fstFile = fullfile(inDir, 'FAST_NTM18B.fst');
tEnd = 660;
outMat = fullfile(fileparts(dirWork), 'dataIn', 'OutDataWind18NTW.mat');
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
    assignin('base', 'OutList', channels'); % as getSimulationOutputTable
    load_system(simMdlname);
    fprintf('FAST run %d of 2 (%s, %g s) ...\n', idxRun, fstFile, tEnd);
    sim(simMdlname, 'SrcWorkspace', 'base');
    close_system(simMdlname, 0);

    % Read the FAST text output (tab delimited, channel names, units)
    [~, fstName] = fileparts(fstFile);
    outFile = fullfile(inDir, [fstName, '.SFunc.out']);
    lines = readlines(outFile);
    idxHead = find(startsWith(strtrim(lines), 'Time'), 1);
    names = strtrim(split(lines(idxHead), char(9)))';
    data = readmatrix(outFile, 'FileType', 'text', 'Delimiter', '\t', ...
        'NumHeaderLines', idxHead + 1);
    data = data(all(isfinite(data), 2), :);
    [isOut, idxCol] = ismember(channels, names);
    assert(all(isOut), 'Channels missing in the FAST output: %s', strjoin(channels(~isOut), ', '));
    OutTableFAST = array2table(data(:, idxCol), 'VariableNames', channels);
    assert(abs(OutTableFAST.Time(end) - tEnd) < DT/2, 'FAST run ended at %g s', OutTableFAST.Time(end));
    windSeed = interp1(OutTableFAST.Time, OutTableFAST.Wind1VelX, time, 'linear', 'extrap');
end

OutTable = OutTableFAST;
OutTable.Time = time; % exact time grid
save(outMat, 'OutTable');
fprintf('Saved %s: %d samples, %g to %g s, wind mean %.2f m/s, std %.2f m/s\n', outMat, ...
    height(OutTable), OutTable.Time(1), OutTable.Time(end), mean(OutTable.Wind1VelX), std(OutTable.Wind1VelX));
