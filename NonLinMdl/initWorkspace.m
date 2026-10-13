%% Clear workspace
clear; close all; clc;
% Remove only this repository's folders from the path (e.g. Sandbox copies).
% restoredefaultpath would also remove user and add-on paths.
dirRepo = fileparts(fileparts(mfilename('fullpath')));
dirPathCell = strsplit(path, pathsep);
dirPathCell = dirPathCell(startsWith(dirPathCell, dirRepo, 'IgnoreCase', ispc));
if ~isempty(dirPathCell)
    rmpath(dirPathCell{:});
end
bdclose all;
clear ABgL_WT; % LUT variables saved as persistent variables
Simulink.data.dictionary.closeAll;

%% Set paths
dataIn = 'dataIn'; % external data from FAST linearization, FASTtool
dataOut = 'dataOut'; % data from simulation runs, saved to not run all simulation again
dirSub = {'subfunctionsMPC', 'subfunctionsUtil', 'WECSCtrl', 'WECSMDl','WECSAct'}; 

% Get relative paths with respect to testCrl.m file
dirFileName = mfilename('fullpath');
dirWork = fileparts(dirFileName);
dirParent = fileparts(dirWork);

% Add parentDir and workdir to path: Can be called in both directories
addpath(dirParent)
addpath(dirWork)

% Add LinMdl (for WECS param in workspace for testing)
addpath(fullfile(dirParent,'LinMdl'));

% Add path to input data (extracted from FAST with FASTtool)
dataInPath = fullfile(dirParent,dataIn);
addpath(dataInPath);

% Add output path to save mat files
dataOutPath = fullfile(dirWork,dataOut); 
if ~isdir(dataOutPath) %#ok<ISDIR>
    mkdir(dataOutPath);
end
addpath(dataOutPath);

% Add subdirerectoreis
for dirIdx = 1: length(dirSub)
    addpath(fullfile(dirWork,dirSub{dirIdx}));
end

%% Load FASTTool wind as test signal
DT = 0.008; % one sample time 
step_thresshold = 0; % per default filter wind
kVeFilt = 0;
kVeFilt_tau = 0.75;
% Rotor-effective wind filter at the qLPV-MPC wind input: time constant
% wecs.Rr*kVeFiltCtrl/V (0: filter off, the measured wind is passed through).
% Same as the filter for the aerodynamic torque in Mdl_*_WindFiltered.
kVeFiltCtrl = 0.75;
% qLPV-MPC weights: stage cost x'Qx + u'Ru, terminal weight P = pMPC*Q
% x = [Tg, beta, OmegaR, dThetaS, dOmegaR, dOmegaG, dTg, dBeta, dV, dTg_ref(k-1), dBeta_ref(k-1)],
% u = [dTg_ref, dBeta_ref]; references [Tg_ss, beta_ss, OmegaR_ss](V_filt)
% Two designs:
% - 'beta' (default): original design (paper/thesis), tracks Tg and beta
% - 'omega': tracks Tg and the rotor speed instead of the pitch (experiments
%   E50, E66, S01-S08). Better on FAST in turbulent wind (the beta design
%   has no rotor speed feedback); needs the large pitch weight, lower values
%   (rMPC(2) <= 1e5) excite the tower. Controller candidate for FAST
%   (experiments/EXPERIMENTS.md, batch U), FAST 18 m/s, 60-600 s: rotor
%   speed std 0.144 rpm (PI 0.286), power std 61 kW (PI 120), tower
%   fore-aft acc. 0.216 m/s^2 (PI 0.190), pitch rate std 1.24 deg/s (PI 0.62).
ctrlDesign = 'beta';
switch ctrlDesign
    case 'omega'
        qMPC = [1 0 1e3 0 1e3 1e3 0 0 0 0 0];
        rMPC = [1 1e6];
        pMPC = 10; % small terminal weight; 10 is more robust over two wind seeds than 1
    otherwise
        qMPC = [1 1e4 0 0 1e3 1e3 0 0 0 0 0];
        rMPC = [1 1e4];
        pMPC = 1e3;
end
% First-order low-pass filters of the measured speeds in the qLPV-MPC [Hz]
% (Inf: no filter): rotor speed for scheduling and the tracked rotor speed,
% rotor and generator speed for the state increments
fcOmegaCtrl = Inf;
fcOmegaStatesCtrl = Inf;
% Sample time of the qLPV-MPC (multiple of DT); the prediction horizon is 40*DTc
DTc = 0.04; % horizon 1.6 s; DTc = DT (0.32 s horizon) is unstable on FAST
load(fullfile(dataInPath,'OutDataSweep.mat'),'OutTable');
%load(fullfile(dataInPath,'OutDataWind18NTM.mat'),'OutTable');

%% Clean-up
clear RESTORE*
clear data*
clear dir*



