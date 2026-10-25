function P = qlmpcDesign(name)
% qlmpcDesign returns the parameters of a stored qLPV-MPC controller design
% (base workspace variables of initWorkspace). The struct can also be passed
% to experiments/run_qlmpc_experiment as params.
%
% Fields: ctrlVariant ('beta': MdlCtrl_qLPVMPC, test_*_qLPVMPCbeta; 'flap':
% MdlCtrl_qLPVMPCflap, test_*_qLPVMPCflap; 'flapIPC': MdlCtrl_qLPVMPCflap
% with IPC, test_FAST_qLPVMPCflapIPC), qMPC, rMPC, pMPC, DTc, fNotchCtrl,
% zetaNotchCtrl, smoothPitchCtrl, towerKFCtrl, kfSigmaVCtrl, kfSigmaACtrl,
% ipcParCtrl ([gainSet, kScale, towerOn] of ipc_dlr; gainSet 0: no IPC),
% kVeFiltCtrl (wind filter at the controller input, 0: off), paperMPCCtrl
% (1: original qLPV-MPC formulation of the paper, HSqLd_WT_beta).
%
% Designs (results: experiments/EXPERIMENTS.md, FAST 18 m/s, 60-600 s,
% rotor speed std / power std / tower fore-aft acc. std / pitch rate std;
% PI: 0.286 rpm / 120 kW / 0.190 m/s^2 / 0.62 deg/s):
% - 'beta': original design (paper/thesis), tracks Tg and beta; good on the
%   Simulink plant models, fails on FAST in turbulent wind
% - 'paper': 'beta' as in the paper (runGeneratePicsPaper): controller
%   sample time DTc = DT = 0.008 s, no wind filter at the controller input,
%   original formulation (Taylor discretization, scheduling)
% - 'omega': tracks Tg and the rotor speed instead of the pitch
%   (0.144 / 61 / 0.216 / 1.24)
% - 'flap': omega design with the tower fore-aft and blade flap prediction
%   model, pitch weight 1e5 (0.106 / 51 / 0.315 / 6.2)
% - 'flapNotch': flap with notch filters at 3P and the drivetrain mode
%   (0.086 / 38 / 0.290 / 3.6)
% - 'flapNotchSmooth': flapNotch with the smooth pitch reference
%   (0.087 / 39 / 0.298 / 3.2)
% - 'flapSmooth12': flap, DTc = 0.12 s, smooth pitch reference
%   (0.150 / 63 / 0.244 / 1.8)
% - 'flapNotchSmooth12': flapSmooth12 with notch filters
%   (0.155 / 66 / 0.233 / 1.65)
% - 'flapNotchSmooth12Tower': flapNotchSmooth12 with the tower velocity
%   weight 1e2 (0.185 / 78 / 0.219 / 1.67)
% - 'flapNotchSmooth12KF', 'flapNotchSmooth12TowerKF': as above with the
%   Kalman filter for the tower and blade flap states (flap_kalman)
% - 'flapNotchIPC': flapNotch with the DLR individual pitch controller on top
%   (ipc_dlr: gain set 2 of initDLRparameter scaled by 0.3, without the
%   tower damping term; ctrlVariant 'flapIPC': test_FAST_qLPVMPCflapIPC,
%   FAST only)
% - 'flapNotchSmooth12IPC': flapNotchSmooth12 with the same IPC on top
%
% qMPC: 'beta' model x = [Tg, beta, OmegaR, dThetaS, dOmegaR, dOmegaG, dTg,
% dBeta, dV, dTg_ref(k-1), dBeta_ref(k-1)]; 'flap' model [Tg beta OmegaR v_t
% | dThetaS dOmegaR dOmegaG dTg dBeta dy_t ddelta dv_t dddelta | dV | dU].

if nargin < 1 || isempty(name)
    name = 'beta';
end

% Defaults: no notch filter, step-wise pitch reference, DTc = 0.04 s
P = struct('ctrlVariant', 'beta', 'qMPC', [], 'rMPC', [], 'pMPC', [], ...
    'DTc', 0.04, 'fNotchCtrl', 0, 'zetaNotchCtrl', 0.1, 'smoothPitchCtrl', 0, ...
    'towerKFCtrl', 0, 'kfSigmaVCtrl', 0.5, 'kfSigmaACtrl', 0.05, 'ipcParCtrl', [0 1 0], ...
    'kVeFiltCtrl', 0.75, 'paperMPCCtrl', 0);
qFlap = [1 0 1e3 0 0 1e3 1e3 0 0 0 0 0 0 0 0 0];
fNotch = [0.605 1.64]; % 3P at rated rotor speed, FAST drivetrain mode

switch name
    case 'beta'
        P.qMPC = [1 1e4 0 0 1e3 1e3 0 0 0 0 0];
        P.rMPC = [1 1e4];
        P.pMPC = 1e3;
    case 'paper'
        P = qlmpcDesign('beta');
        P.DTc = 0.008;
        P.kVeFiltCtrl = 0;
        P.paperMPCCtrl = 1;
    case 'omega'
        P.qMPC = [1 0 1e3 0 1e3 1e3 0 0 0 0 0];
        P.rMPC = [1 1e6];
        P.pMPC = 10;
    case {'flap', 'flapNotch', 'flapNotchSmooth', 'flapSmooth12', 'flapNotchSmooth12', ...
            'flapNotchSmooth12Tower', 'flapNotchSmooth12KF', 'flapNotchSmooth12TowerKF'}
        P.ctrlVariant = 'flap';
        P.qMPC = qFlap;
        P.rMPC = [1 1e5];
        P.pMPC = 10;
        if contains(name, 'Notch'), P.fNotchCtrl = fNotch; end
        if contains(name, 'Smooth'), P.smoothPitchCtrl = 1; end
        if contains(name, '12'), P.DTc = 0.12; end
        if contains(name, 'Tower'), P.qMPC(4) = 1e2; end % tower velocity v_t
        if endsWith(name, 'KF'), P.towerKFCtrl = 1; end
    case {'flapNotchIPC', 'flapNotchSmooth12IPC'}
        P = qlmpcDesign(erase(name, 'IPC'));
        P.ctrlVariant = 'flapIPC';
        P.ipcParCtrl = [2 0.3 0];
    otherwise
        error('Unknown qLPV-MPC design %s.', name);
end
end
