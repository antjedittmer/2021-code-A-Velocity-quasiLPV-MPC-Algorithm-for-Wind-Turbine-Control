# qLPV-MPC experiments on FAST (rotor speed tracking)

Goal: make the qLPV-MPC work on FAST v8.16 (NREL 5MW, 18 m/s turbulent wind) at
least as well as the PI baseline, in particular by tracking the rotor speed.

## Setup

- `run_qlmpc_experiment.m`: one closed-loop run with parameter overrides of the
  base workspace (see `initWorkspace.m`). Result: `results/<name>.mat`
  (metrics `M`, parameters, time series `T5` with 0.04 s resolution), one line
  in `results_log.txt`.
- `run_queue.py`: runs the experiments of a queue file (`queue_batch*.txt`) in
  parallel MATLAB processes; logs in `logs/`.
- Plants: `FAST` (FAST v8.16, `FASTinputfiles/FAST_NTW18.fst`, TurbSim NTM class A,
  hub height mean 18.2 m/s, TI 16 %), `SimFilt` (`Mdl_BianchiOL_WindFiltered`),
  `SimPlain` (`Mdl_BianchiOL`).
- Metrics for t >= 15 s; status `failed` = simulation aborted (FAST diverged),
  `unstable` = rotor speed above 15 rpm.

PI baseline on FAST (same wind, `dataOut/OutTableTest2FASTNTW18.mat`):
rotor speed 11.98 ± 0.29 rpm (max 12.89), power 4999 ± 120 kW,
tower fore-aft acc. 0.191 m/s², pitch rate std 0.63 deg/s.

## Controller parameters (base workspace, `initWorkspace.m`)

| Variable | Meaning | Default |
|---|---|---|
| `qMPC` | state weights `[Tg, beta, OmegaR, dThetaS, dOmegaR, dOmegaG, dTg, dBeta, dV, dU(k-1)]` | `[1 1e4 1e3 0 1e3 1e3 0 0 0 0 0]` |
| `rMPC` | input weights `[dTg, dBeta]` | `[1 1e4]` |
| `pMPC` | terminal weight factor (P = pMPC*Q) | `1e3` |
| `kVeFiltCtrl` | rotor-effective wind filter, tau = 63 m * kVeFiltCtrl / V | `0.75` |
| `fcOmegaCtrl` | low-pass [Hz] of the rotor speed for scheduling and the tracked speed | `Inf` (off) |
| `fcOmegaStatesCtrl` | low-pass [Hz] of rotor and generator speed for the state increments | `Inf` (off) |
| `DTc` | controller sample time (horizon 40*DTc) | `DT` = 0.008 s |

## Batch 1: weights at DTc = 8 ms (horizon 0.32 s)

| Name | Change | Result on FAST |
|---|---|---|
| E01 | old controller (omega weight 0) | runs, 11.44 ± 1.26 rpm, max 15.1 |
| E02 | beta 1e4, omega 1e3 | runs, 11.84 ± 0.44 rpm, max 13.5, pitch rate std 9.5 deg/s |
| E09 | as E02, speed filter 0.25 Hz | runs, 11.83 ± 0.45 rpm, tower -14 %, pitch rate std 8.2 deg/s |
| E03-E08, E10 | beta weight 0 (omega 1e2..1e4, with/without filter, rMPC(2) = 1e6) | all diverge (FAST aborts after 15-108 s) |

With a horizon of only 0.32 s, tracking the rotor speed without the beta
reference is unstable.

## Bug: discretization of the prediction model

`ABqLd_WT_beta` discretized the model with a 4th order Taylor series. For the
torque actuator (kappa_g = 0.01 s) this gives exp(-DT/kappa) = 5.0 instead of
0.018 at DT = 0.04 s and 110 at 0.08 s (at 8 ms: 0.452 instead of 0.449).
Replaced by the exact zero-order-hold discretization (`expm`). Batch 2 (DTc =
0.04/0.08 s with the Taylor series) is invalid and was stopped.

## Bug: invalid first samples at DTc > DT

With DTc > DT, Simulink inserts rate transitions whose first output is 0. With
V = 0 or OmegaR = 0 the cQ table lookup gives NaN, the NaN stays in the
predicted trajectory Xk, and the controller drives pitch and torque to their
lower limits at the maximum rate (min(max(NaN, -lim), lim) = -lim). Fixed by
guards in the controller: scheduling values limited to plausible ranges,
non-finite trajectories restarted from the current state, non-finite
solutions treated as solver failure (and limited in HSqLd_WT_beta).
Regression: E11 (after the DTc change, DTc = DT) is identical to E02; E21
(expm) differs only marginally from E02.

## Batches 3-5: longer horizon (DTc = 0.04 s, horizon 1.6 s)

| Name | Change | RotSpd [rpm] | Power std [kW] | Tower FA [m/s^2] | Pitch rate std [deg/s] |
|---|---|---|---|---|---|
| PI | | 11.98 ± 0.29, max 12.89 | 120 | 0.191 | 0.63 |
| E31 | beta 1e4, omega 1e3 | 12.06 ± 0.28, max 13.09 | 128 | 0.458 | 10.1 |
| E35 | as E31, rMPC(2) 1e6 | 12.06 ± 0.15, max 12.76 | 75 | 0.739 | 7.8 |
| E32/E33/E38 | beta 0, rMPC(2) 1e4 (DTc 0.04/0.08) | unstable | | | |
| E37 | as E31, DTc 0.08 | 11.86 ± 0.42, max 13.15 | 177 | 0.390 | 9.0 |
| E39 | beta 0, omega 1e3, rMPC(2) 1e6 | 12.01 ± 0.12, max 12.47 | 54 | 0.489 | 3.9 |
| E44 | E39, rMPC(2) 1e7 | 12.00 ± 0.10, max 12.35 | 43 | 0.252 | 2.1 |
| E47 | E39, fcOmegaCtrl 0.25 Hz | 12.00 ± 0.11, max 12.45 | 48 | 0.272 | 3.2 |
| E50 | E39, pMPC 1 (instead of 1000) | 11.98 ± 0.18, max 12.55 | 77 | 0.207 | 1.0 |
| E41 | E31, also filter the state increments | 11.98 ± 0.29, max 13.04 | 122 | 1.441 | 7.5 |

Findings:
- Replacing beta by omega (beta weight 0) works on FAST with the longer
  horizon and a strong penalty on the pitch increments (rMPC(2) >= 1e6).
- The pitch "chatter" consists of 3P (about 0.6 Hz), the drivetrain mode
  (1.64 Hz) and, with tight speed tracking, the first tower fore-aft mode
  (0.32 Hz): pitching at the tower frequency reduces the tower damping
  (E39: tower acceleration in 0.25-0.4 Hz 15 times the PI value). The
  prediction model contains no tower.
- The terminal weight pMPC = 1000 makes the controller aggressive; pMPC = 1
  (E50) brings the tower acceleration almost to the PI level.
- Filtering the speed increments of the state (fcOmegaStatesCtrl) is harmful.
- The computer was in sleep mode from 00:22 to 03:38 during batch 5 (runs
  paused, results not affected).
