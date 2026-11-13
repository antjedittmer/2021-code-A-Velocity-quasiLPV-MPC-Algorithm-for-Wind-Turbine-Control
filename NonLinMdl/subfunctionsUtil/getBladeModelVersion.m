function [bladeVersion, bladeCompare] = getBladeModelVersion()
% getBladeModelVersion returns the version of the blade model of Mdl2 set
% by the base workspace variable bladeModelVersion, and whether the
% original Mdl2 is shown additionally in the plots (base workspace variable
% bladeModelCompare, review of the blade model).
%
% - 'original' (default, roll-back): rigid blade with modal mass, stiffness
%   and thrust load point at r_b = 0.75 R (initModel5MWNREL)
% - 'B1': variant B, modal mass at the fixed reference radius r_ref, thrust
%   load point r_b(lambda, beta) and stiffness Kb(lambda, beta) fitted to
%   the quasi-steady FAST blade tip deflection (fitBladeModelB)
% - 'B2': as B1, modal mass scaled with Kb (constant blade frequency)

bladeVersion = 'original';
bladeCompare = false;
try %#ok<TRYNC>
    if evalin('base', 'exist(''bladeModelVersion'', ''var'')')
        bladeVersion = evalin('base', 'bladeModelVersion');
    end
    if evalin('base', 'exist(''bladeModelCompare'', ''var'')')
        bladeCompare = logical(evalin('base', 'bladeModelCompare'));
    end
end
if ~any(strcmp(bladeVersion, {'original', 'B1', 'B2'}))
    error('Unknown blade model version %s (original, B1 or B2).', bladeVersion);
end
bladeCompare = bladeCompare && ~strcmp(bladeVersion, 'original');
end
