function [modelVersion, stepDataMat, dataOutSub] = getWecsModelVersion()
% getWecsModelVersion returns the version of the simplified WECS models
% (Mdl1, Mdl2) set by the base workspace variable wecsModelVersion
% (initWorkspace), the corresponding FAST step reference and the subfolder
% of NonLinMdl/dataOut for the simulation results ('' or 'paper').
%
% - 'FAST' (default): models matched to the FAST v8.16 S-function of this
%   repository: effective wind factor wecs.kWindEff = 0.985, local +5 % Cq
%   correction near rated, step reference OutDataStep.mat (FAST v8.16)
% - 'paper': original models of the paper/thesis (runGeneratePicsPaper):
%   no effective wind factor (kWindEff = 1), cos(5 deg) in the wind filter
%   (Lib_WindFilter), no Cq correction, step reference
%   OutDataStepFASTTool.mat (FASTTool)
%
% The data dictionaries are updated with updateDDMdl1.

modelVersion = 'FAST';
try %#ok<TRYNC>
    if evalin('base', 'exist(''wecsModelVersion'', ''var'')')
        modelVersion = evalin('base', 'wecsModelVersion');
    end
end

switch modelVersion
    case 'FAST'
        stepDataMat = 'OutDataStep.mat';
        dataOutSub = '';
    case 'paper'
        stepDataMat = 'OutDataStepFASTTool.mat';
        dataOutSub = 'paper';
    otherwise
        error('Unknown WECS model version %s (FAST or paper).', modelVersion);
end
end
