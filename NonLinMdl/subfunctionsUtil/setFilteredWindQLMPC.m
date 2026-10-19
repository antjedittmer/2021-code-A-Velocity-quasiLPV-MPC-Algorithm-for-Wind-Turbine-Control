function oldValues = setFilteredWindQLMPC(simMdlname)
% setFilteredWindQLMPC sets the rotor-effective wind filter for a
% qLPV-MPC test model (turbulent wind):
% - plant: the WECS model reference of simMdlname is set to
%   Mdl_BianchiOL_WindFiltered (in memory only, close the model without
%   saving) with the filter settings of runCompareModels (step_thresshold,
%   kVeFilt, kVeFilt_tau in the base workspace)
% - controller: the qLPV-MPC is scheduled with the same filtered wind as the
%   aerodynamic torque of the plant, i.e. the wind is projected onto the
%   rotor axis (kWindEffCtrl = wecs.kWindEffFilt, cos(5 deg) in the paper
%   version) and filtered with kVeFiltCtrl = kVeFilt_tau
%
% oldValues: base workspace values before the call (struct), restore them
% with restoreBaseVars(oldValues).

% Projection factor of the plant wind filter (wecs in DD_Mdl1.sldd)
ddObj = Simulink.data.dictionary.open('DD_Mdl1.sldd');
wecs = getValue(getEntry(getSection(ddObj, 'Design Data'), 'wecs'));

newValues = struct('step_thresshold', 1, 'kVeFilt', 0.25, 'kVeFilt_tau', 0.75, ...
    'kVeFiltCtrl', 0.75, 'kWindEffCtrl', wecs.kWindEffFilt);

names = fieldnames(newValues);
oldValues = struct;
for idx = 1:numel(names)
    oldValues.(names{idx}) = evalin('base', names{idx});
    assignin('base', names{idx}, newValues.(names{idx}));
end

load_system(simMdlname);
blk = find_system(simMdlname, 'BlockType', 'ModelReference', 'Name', 'WECS Model');
set_param(blk{1}, 'ModelName', 'Mdl_BianchiOL_WindFiltered');
end
