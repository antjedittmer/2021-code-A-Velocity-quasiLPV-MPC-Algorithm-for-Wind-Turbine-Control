function [timeBlade, thetaB] = getBladeStateMdl2(matFileBlade, loadData, OutTable, simMdlname2, OutTableTest2)
% getBladeStateMdl2 returns the blade state theta_b (rad) of Mdl2
% (Mdl_BianchiOL, second generalized coordinate q(2) of the Lagrange
% model: rigid blade rotation, out-of-plane displacement r*theta_b at
% radius r) for the plots of runCompareModels.
%
% The blade state is not part of the stored simulation output tables
% (OutTableTest2*.mat). It is loaded from matFileBlade if available (and
% loadData), otherwise Mdl2 is simulated again with the states logged
% (SaveState, no change of the Simulink models) and the result is saved in
% matFileBlade. The stored output table OutTableTest2 is not changed; the
% outputs of the new simulation are compared with it (maximum relative
% difference in the console).
%
% Inputs
% - matFileBlade: mat file of the blade state
% - loadData: load the blade state if available
% - OutTable: FAST simulation table (wind input of the Simulink model)
% - simMdlname2: Simulink test model of Mdl2 (current WECS model reference)
% - OutTableTest2: stored output table of Mdl2 (for the comparison)
%
% Outputs
% - timeBlade: time vector (s)
% - thetaB: blade state theta_b (rad)

if exist(matFileBlade, 'file') == 2 && loadData
    tmp = load(matFileBlade, 'timeBlade', 'thetaB');
    timeBlade = tmp.timeBlade;
    thetaB = tmp.thetaB;
    return
end

WS = 'base';
assignin(WS, 'OutTable', OutTable)
load_system(simMdlname2);
% LimitDataPoints off: the test models keep only the last 1000 points
out = sim(simMdlname2, 'SrcWorkspace', WS, 'SaveState', 'on', ...
    'StateSaveName', 'xoutBlade', 'SaveFormat', 'Dataset', 'LimitDataPoints', 'off');

% Integrator1 of the WECS subsystem: generalized coordinates
% q = [x_t, theta_b, y_t, phi_r, phi_g]
xout = out.xoutBlade;
idxQ = [];
for idx = 1:xout.numElements
    blkPath = xout{idx}.BlockPath.getBlock(xout{idx}.BlockPath.getLength);
    if endsWith(blkPath, '/WECS/Integrator1')
        idxQ = idx;
        break
    end
end
if isempty(idxQ)
    error('getBladeStateMdl2:noState', 'State WECS/Integrator1 not found in %s.', simMdlname2);
end
q = xout{idxQ}.Values;
timeBlade = q.Time(:);
qData = squeeze(q.Data);
if size(qData, 1) ~= numel(timeBlade)
    qData = qData.';
end
thetaB = qData(:, 2);
if timeBlade(1) > OutTable.Time(1) + 1 || timeBlade(end) < OutTable.Time(end) - 1
    error('getBladeStateMdl2:timeRange', ...
        'Blade state logged only from %.1f s to %.1f s.', timeBlade(1), timeBlade(end));
end

% Comparison with the stored outputs (same simulation, outputs unchanged)
nOut = min(height(OutTableTest2), size(out.OutDataTest, 1));
nVar = min(width(OutTableTest2), size(out.OutDataTest, 2));
stored = OutTableTest2{1:nOut, 1:nVar};
newOut = out.OutDataTest(1:nOut, 1:nVar);
relDiff = max(abs(newOut - stored), [], 1) ./ max(max(abs(stored), [], 1), eps);
fprintf('%s: maximum relative difference of the new Mdl2 run to the stored outputs: %.2e\n', ...
    matFileBlade, max(relDiff));

save(matFileBlade, 'timeBlade', 'thetaB');
end
