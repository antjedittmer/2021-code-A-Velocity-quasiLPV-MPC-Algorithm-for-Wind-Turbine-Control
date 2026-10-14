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
% Stored designs (subfunctionsMPC/qlmpcDesign.m, results in
% experiments/EXPERIMENTS.md): 'beta' (default, original paper/thesis
% design), 'omega' (rotor speed tracking, candidate for FAST), 'flap',
% 'flapNotch', 'flapNotchSmooth', 'flapSmooth12', 'flapNotchSmooth12'
% (tower fore-aft and blade flap prediction model; need the controller
% MdlCtrl_qLPVMPCflap, i.e. test_SimulinkMdl2_qLPVMPCflap or
% test_FAST_qLPVMPCflap). ctrlVariant names the controller model variant.
ctrlDesign = 'beta';
designMPC = qlmpcDesign(ctrlDesign);
ctrlVariant = designMPC.ctrlVariant;
qMPC = designMPC.qMPC;
rMPC = designMPC.rMPC;
pMPC = designMPC.pMPC;
% Sample time of the qLPV-MPC (multiple of DT); the prediction horizon is 40*DTc
DTc = designMPC.DTc; % 0.04: horizon 1.6 s; DTc = DT (0.32 s horizon) is unstable on FAST
% Only MdlCtrl_qLPVMPCflap: notch filters in the speed filters at the
% frequencies fNotchCtrl [Hz] (0: off; e.g. [0.605 1.64] for 3P at rated
% speed and the FAST drivetrain mode) with the pole damping zetaNotchCtrl,
% and the pitch reference ramped at DT between the controller updates
% (smoothPitchCtrl = 1) instead of step-wise (0)
fNotchCtrl = designMPC.fNotchCtrl;
zetaNotchCtrl = designMPC.zetaNotchCtrl;
smoothPitchCtrl = designMPC.smoothPitchCtrl;
% Only MdlCtrl_qLPVMPCflap: tower and blade flap states from a Kalman filter
% (towerKFCtrl = 1, flap_kalman; process noise: hub/rotor-effective wind
% difference kfSigmaVCtrl [m/s], measurement noise of the tower acceleration
% kfSigmaACtrl [m/s^2]) instead of the leaky integration of the tower
% acceleration (corner frequency fcTowerCtrl [Hz]) and flap_observer (0)
towerKFCtrl = designMPC.towerKFCtrl;
kfSigmaVCtrl = designMPC.kfSigmaVCtrl;
kfSigmaACtrl = designMPC.kfSigmaACtrl;
fcTowerCtrl = 0.02;
towerSignCtrl = 1; % sign of the measured tower acceleration in the controller
% First-order low-pass filters of the measured speeds in the qLPV-MPC [Hz]
% (Inf: no filter): rotor speed for scheduling and the tracked rotor speed,
% rotor and generator speed for the state increments
fcOmegaCtrl = Inf;
fcOmegaStatesCtrl = Inf;
load(fullfile(dataInPath,'OutDataSweep.mat'),'OutTable');
%load(fullfile(dataInPath,'OutDataWind18NTM.mat'),'OutTable');

%% Clean-up
clear RESTORE*
clear data*
clear dir*



