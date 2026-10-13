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

## Batch 6: ratio of speed weight and pitch increment weight

With beta weight 0 the result depends mainly on the ratio q_omega / r_beta
(E54, E55, E56 have the same ratio 3e-5 and practically identical results).
FAST, DTc = 0.04 s, pMPC = 1000:

| q_omega / r_beta | RotSpd std / max [rpm] | Power std [kW] | Tower FA [m/s^2] | Pitch rate std [deg/s] |
|---|---|---|---|---|
| 1e-3 (E39) | 0.12 / 12.47 | 54 | 0.489 | 3.9 |
| 1e-4 (E44) | 0.10 / 12.35 | 43 | 0.252 | 2.1 |
| 3e-5 (E52, E54-E56) | 0.14-0.15 / 12.43 | 60-63 | 0.215 | 1.3 |
| 1e-5 (E53) | 0.29 / 12.88 | 120 | 0.194 | 0.73 |
| PI | 0.29 / 12.89 | 120 | 0.191 | 0.63 |

For small ratios the qLPV-MPC approaches the PI behaviour; pMPC = 1 (E50)
gives a better trade-off.

## Tower model (Mdl2 without blade flap and tower side-side)

`MdlCtrl_qLPVMPCtower` (test models `test_FAST_qLPVMPCtower`,
`test_SimulinkMdl2_qLPVMPCtower`, runner option `'ctrlVariant','tower'`):
prediction model `ABqL_WT_tower` = Mdl2 (`initModel5MWNREL`) without the blade
flap and the tower side-side states,
x = [ThetaS, OmegaR, OmegaG, Tg, beta, x_t, v_t],
mtb*dv_t = -Kt*x_t - Bt*v_t + Ft(V - v_t, OmegaR, beta), aerodynamic torque
with the relative wind V - v_t, thrust sensitivities from the cT table.
Outputs y = [Tg, beta, OmegaR, v_t] (v_t with reference 0 for tower damping);
qMPC has 14 elements
[Tg beta OmegaR v_t | dThetaS dOmegaR dOmegaG dTg dBeta dx_t dv_t | dV | dU].
The tower velocity and displacement increments are estimated from the
measured fore-aft acceleration by leaky integration (`fcTowerCtrl` = 0.02 Hz).
At 18 m/s the model has the tower mode at 0.318 Hz (damping 5.2 %) and the
drivetrain mode at 2.22 Hz.

## Batch 7: small terminal weight (pMPC = 1)

FAST, DTc = 0.04 s, beta weight 0:

| Name | Change | RotSpd [rpm] | Power std [kW] | Tower FA [m/s^2] | Pitch rate std [deg/s] |
|---|---|---|---|---|---|
| E50 / E66 | omega 1e3, rMPC(2) 1e6, pMPC 1 / 0 | 11.98 ± 0.18-0.19, max 12.55 | 77-80 | 0.205 | 1.0 |
| E64 | pMPC 10 | 11.98 ± 0.14, max 12.44 | 61 | 0.218 | 1.2 |
| E65 | omega 3e3 | 11.98 ± 0.10, max 12.34 | 43 | 0.257 | 1.9 |
| E63 | speed filter 0.25 Hz | 11.97 ± 0.27, max 12.88 | 112 | 0.198 | 1.05 |
| E62 | rMPC(2) 1e7 | 11.87 ± 0.83, max 14.53 | 345 | 0.182 | 0.36 |
| E67 / E68 | rMPC(2) 1e5 / omega 1e4 | 12.00 ± 0.32 | 134 | 1.40-1.44 | 7.1 |
| E69 | DTc 0.08 s | 11.97 ± 0.21, max 12.63 | 87 | 0.204 | 1.0 |
| E71 | beta weight 1e4 added | 11.57 ± 0.98, max 14.60 | 408 | 0.201 | 0.79 |

All variants without a tower model lie on one trade-off curve between speed
regulation and tower load; the weights select the point on the curve. The
beta reference (E71) brings back the speed offset: beta weight 0 is the
right choice.

Wind sweep on FAST (E60: beta 0, omega 1e3, rMPC(2) 1e7, pMPC 1000) against
the PI: below rated practically equal; above rated (13-23 m/s) rotor speed
11.99 ± 0.05 rpm (PI 12.03 ± 0.11), max 12.19 (PI 12.32), tower FA 0.086
(PI 0.061) m/s^2, pitch rate std 0.45 (PI 0.23) deg/s.

## Batches 8-10: tower model

| Name | Plant | Tower weight v_t | Result |
|---|---|---|---|
| T01 | FAST | 0 | 11.98 ± 0.19 rpm, tower 0.206, pitch rate 0.99 (as E50) |
| T10 / T11 / T12 | FAST | 1e1 / 1e2 / 3e2 | tower 0.207 / 0.228 / 0.549 |
| T02-T06 | FAST | 1e3 ... 1e6 | diverge (tower oscillation at 0.29 Hz) |
| T00 / T13 / T24 | SimFilt | 0 / 1e2 / 1e3 | tower 0.139 / 0.131 / 0.151 (speed std 0.22 / 0.26 / 0.48) |
| T20-T22 | FAST | sign of the acceleration inverted | worse (1e3, 1e4 diverge; 1e2: tower 0.273) |
| T23 | SimFilt | inverted sign, 1e3 | tower 0.190 (positive sign: 0.151) |

The tower model works (T01 = E50), and on the Simulink plant (Mdl2) the
tower velocity weight damps the tower. On FAST it increases the tower
acceleration with either sign of the measurement: the reduced tower model
(without blade flexibility) does not represent the FAST tower/rotor coupling
well enough for active damping. Not pursued further.

## Recommendation (new defaults in initWorkspace.m)

beta weight 0 (rotor speed tracking instead of the beta reference),
`qMPC = [1 0 1e3 0 1e3 1e3 0 0 0 0 0]`, `rMPC = [1 1e6]`, `pMPC = 1`,
`DTc = 0.04` s, `kVeFiltCtrl = 0.75`, no speed filters.
On FAST (18 m/s, E50): rotor speed 11.98 ± 0.18 rpm (PI ± 0.29), max 12.55
(PI 12.89), power std 77 kW (PI 120), tower fore-aft acceleration 0.207 m/s^2
(PI 0.191), pitch rate std 1.0 deg/s (PI 0.63). More speed precision is
available along the trade-off curve (e.g. E65: ± 0.10 rpm, tower 0.257).

Not tried yet: a notch filter for 3P / drivetrain mode in the speed
measurement, a better tower/blade model, estimating the rotor-effective
wind speed instead of measuring the hub wind.

## Second wind field (seed B)

FAST, `FASTinputfiles/FAST_NTW18b.fst` (TurbSim NTM class A, new seed,
`wind18_TurbKlasseA_NewSeed.bts`, hub height mean 17.8 m/s, TI 15.4 %; the
.bts file is not in the repository). PI: `dataOut/OutTableTest2FASTNTW18b.mat`.
S04 (defaults on seed A) reproduces E50 exactly.

| | RotSpd [rpm] | Power std [kW] | Tower FA [m/s^2] | Pitch rate std [deg/s] |
|---|---|---|---|---|
| PI | 11.98 ± 0.25, max 12.69 | 105 | 0.200 | 0.61 |
| S01 defaults (E50 setting) | 11.98 ± 0.20, max 13.23 | 86 | 0.215 | 1.0 |
| S02 E65 setting (omega 3e3) | 11.98 ± 0.11, max 12.76 | 47 | 0.278 | 1.9 |
| S03 E44 setting (rMPC(2) 1e7, pMPC 1e3) | 12.00 ± 0.10, max 12.49 | 43 | 0.280 | 2.1 |

The default setting generalizes only partly (speed std 19 % below the PI,
but a higher maximum speed). The faster settings E44/E65 reduce the speed and
power scatter by a factor of about 2.5 on both seeds, at about 40 % more
tower fore-aft acceleration and three times the pitch activity. The pitch
rate excess comes mainly from 2-6 Hz (step-wise reference at DTc = 0.04 s),
3P and the drivetrain mode.

## Final tuning over both seeds (S05-S08)

| | Seed A: RotSpd std / max [rpm], tower [m/s^2] | Seed B: RotSpd std / max [rpm], tower [m/s^2] | Pitch rate std [deg/s] |
|---|---|---|---|
| PI | 0.29 / 12.89, 0.191 | 0.25 / 12.69, 0.200 | 0.6 |
| pMPC 1 (E50, S01) | 0.18 / 12.55, 0.207 | 0.20 / 13.23, 0.215 | 1.0 |
| **pMPC 10 (S05, S06)** | **0.14 / 12.44, 0.218** | **0.15 / 12.91, 0.230** | 1.2 |
| omega 2e3, pMPC 1 (S07, S08) | 0.12 / 12.38, 0.234 | 0.13 / 12.88, 0.246 | 1.5 |

New default: pMPC = 10 (speed and power scatter about half of the PI on
both seeds, maximum speed about as the PI, tower fore-aft acceleration
+14-15 %, pitch activity twice the PI).
Note: dataOut/OutTableMPC*.mat (runCompareCtrl) were generated with pMPC = 1.

## Evaluation window

From now on the metrics are evaluated over 15-400 s for the Simulink plant
models (fast decaying transients) and over 60-600 s for FAST (runCompareCtrl,
run_qlmpc_experiment, summarize_experiments). The tables above were computed
for t >= 15 s; results_summary.md re-evaluates all experiments from the stored
time series (0.04 s resolution) with the new windows. The rankings do not
change, e.g. FAST 18 m/s, 60-600 s:

| | Seed A: RotSpd std / max [rpm], tower [m/s^2] | Seed B: RotSpd std / max [rpm], tower [m/s^2] | Pitch rate std [deg/s] |
|---|---|---|---|
| PI | 0.286 / 12.89, 0.190 | 0.246 / 12.68, 0.200 | 0.62 |
| default (pMPC 10, S05/S06) | 0.144 / 12.44, 0.215 | 0.155 / 12.91, 0.229 | 1.17 |
| pMPC 1 (E50/S01) | 0.182 / 12.55, 0.205 | 0.208 / 13.23, 0.216 | 0.96 |

runCompareCtrl (std ratio PI/qLPV-MPC, > 1: qLPV-MPC better), FAST 18 m/s
60-600 s: tower 0.88, power tracking error 1.96, pitch rate 0.50; FAST sweep
60-400 s: tower 0.91, power 1.09, pitch rate 0.63.

## Batch U: tower model with the omega design after kWindEff and the Cq correction

Omega weights, pMPC = 10, tower variant `qMPC = [1 0 1e3 vt 0 1e3 1e3 0 ...]`
(14 elements), 18 m/s, metrics as above (SimFilt 15-400 s, FAST 60-600 s).
SimFilt PI (dataOut/OutTableTest2SimFiltNTW18.mat): speed std 0.375 rpm,
power std 157 kW, tower 0.125 m/s^2, pitch rate 0.72 deg/s.

| Name | Plant | vt | rMPC(2) | DTc | RotSpd std / max | Power std | Tower FA | Pitch rate std |
|---|---|---|---|---|---|---|---|---|
| U00 (no tower model) | SimFilt | - | 1e6 | 0.04 | 0.186 / 12.51 | 78 | 0.141 | 0.99 |
| U vt0 | SimFilt | 0 | 1e6 / 1e5 | 0.04 | 0.189 / 0.476 | 79 / 201 | 0.140 / 0.979 | 0.98 / 7.8 |
| U vt1e2 | SimFilt | 1e2 | 1e6 / 1e5 | 0.04 | 0.210 / 1.196 | 88 / 500 | 0.134 / 0.968 | 0.91 / 6.1 |
| U vt1e3 | SimFilt | 1e3 | 1e6 / 1e5 | 0.04 | 0.349 / 0.681 | 146 / 285 | 0.165 / 0.409 | 0.79 / 2.3 |
| DTc08 | SimFilt | 0 / 1e2 | 1e5 | 0.08 | 0.144 / 0.127 | 61 / 53 | 0.214 / 0.208 | 1.76 / 1.49 |
| DTc12 | SimFilt | 0 / 1e2 / 1e3 / 3e3 | 1e5 | 0.12 | 0.150-0.155 | 63-65 | 0.152 / 0.147 / 0.170 / 0.197 | 1.32 / 1.15 / 1.01 / 0.98 |
| A (SimFilt) | SimFilt | 1e2 | 1e6 | 0.12 | 0.297 / 12.78 | 124 | 0.120 | 0.74 |
| C (SimFilt) | SimFilt | 1e3 | 1e5 | 0.16 | 0.188 / 12.52 | 79 | 0.141 | 0.92 |
| rb3e4 | SimFilt | 1e2 / 1e3 | 3e4 | 0.12 / 0.16 | 0.132 / 0.136 | 56 / 57 | 0.231 / 0.173 | 1.90 / 1.13 |
| U_FAST_omega_noTower | FAST | - | 1e6 | 0.04 | 0.144 / 12.44 | 61 | 0.216 | 1.24 |
| U_FAST_A | FAST | 1e2 | 1e6 | 0.12 | 0.255 / 12.72 | 107 | 0.204 | 1.07 |
| U_FAST_B | FAST | 1e2 | 1e5 | 0.12 | 0.341 / 14.35 | 142 | 0.409 | 5.5 |
| U_FAST_C | FAST | 1e3 | 1e5 | 0.16 | 0.472 / 14.47 | 198 | 0.613 | 3.5 |

- With the 1.6 s horizon (DTc = 0.04) a lower pitch weight (1e5) gives a
  limit cycle at 0.24 Hz (rotor speed loop coupled with the tower, below the
  tower mode at 0.34 Hz). A longer horizon (DTc = 0.12, 4.8 s) removes it on
  the Simulink plant.
- On FAST the low pitch weight still destabilizes the tower mode (0.32 Hz)
  even with the tower model and the long horizon (B, C). The reduced tower
  model does not represent the FAST tower/rotor coupling well enough.
- Tower weights above 1e2 do not reduce the tower acceleration.
- Best FAST trade-off of this batch: A (between the PI and the omega design
  without tower model in all metrics).

## Controller candidates

| Name | Selection | Use |
|---|---|---|
| beta (default) | `ctrlDesign = 'beta'` in initWorkspace: `qMPC = [1 1e4 0 0 1e3 1e3 0 0 0 0 0]`, `rMPC = [1 1e4]`, `pMPC = 1e3` | original paper/thesis design; Simulink plant models (fails on FAST in turbulent wind) |
| omega (candidate for FAST) | `ctrlDesign = 'omega'`: `qMPC = [1 0 1e3 0 1e3 1e3 0 0 0 0 0]`, `rMPC = [1 1e6]`, `pMPC = 10`, `DTc = 0.04` | FAST 18 m/s (U_FAST_omega_noTower): speed std 0.144 rpm, power std 61 kW, tower 0.216 m/s^2, pitch rate 1.24 deg/s (PI 0.286 / 120 / 0.190 / 0.62) |
| flap (candidate) | `ctrlDesign = 'flap'`: MdlCtrl_qLPVMPCflap (test_*_qLPVMPCflap), `qMPC = [1 0 1e3 0 0 1e3 1e3 0 0 0 0 0 0 0 0 0]`, `rMPC = [1 1e5]`, `pMPC = 10`, `DTc = 0.04` | FAST 18 m/s (F_FAST_vt0_rb1e5_DTc04): speed std 0.106 rpm, power std 51 kW, tower 0.315 m/s^2, pitch rate 6.2 deg/s (drivetrain mode and 3P); SimFilt: 0.077 / 33 / 0.178 / 1.12 |

## Model check at the tower mode (check_tower_fit.m)

Frequency response pitch -> tower fore-aft acceleration and rotor speed at
12/18/22 m/s (linearized FAST vs. Mdl1, Mdl2, the tower model of the
qLPV-MPC and the new model with blade flap). The magnitudes agree within
about +-20 % around the tower peak (as in the thesis Bode plots, which show
only magnitudes), but the phase of the tower model (and Mdl1) leads FAST by
15-27 deg (acceleration) and 18-35 deg (rotor speed) in 0.2-0.32 Hz, the band
of the FAST oscillations. Mdl2 and the blade flap model ABqL_WT_flap stay
within -8...+8 deg (acceleration) and -6...-12 deg (rotor speed).

## Batch F: prediction model with blade flap (MdlCtrl_qLPVMPCflap)

`ABqL_WT_flap` = Mdl2 without tower side-side, x = [ThetaS, OmegaR, OmegaG,
Tg, beta, y_t, delta, v_t, ddelta]; velocity form `ABqLd_WT_flap`,
`HSqLd_WT_flap` (16 augmented states, qMPC with 16 elements
[Tg beta OmegaR v_t | dThetaS dOmegaR dOmegaG dTg dBeta dy_t ddelta dv_t dddelta | dV | dU]).
The blade flap state is not measured: `flap_observer` integrates the Mdl2
blade equation driven by the measured tower acceleration and velocity and the
deviations of rotor speed, pitch and wind from their 0.02 Hz means (stable
through the aerodynamic blade damping, about 66 % at 18 m/s; correlation
with the true state in the linear model 0.9-1.0). Runner option
`'ctrlVariant','flap'`. Omega weights, pMPC = 10, 18 m/s:

| Name | Plant | vt | rMPC(2) | DTc | RotSpd std / max | Power std | Tower FA | Pitch rate std |
|---|---|---|---|---|---|---|---|---|
| F vt0 rb1e6 | SimFilt | 0 | 1e6 | 0.04 | 0.212 / 12.60 | 89 | 0.122 | 0.75 |
| F vt0 rb1e5 | SimFilt | 0 | 1e5 | 0.04 | 0.077 / 12.26 | 33 | 0.178 | 1.12 |
| F vt0 rb1e5 | SimFilt | 0 | 1e5 | 0.12 | 0.140 / 12.38 | 59 | 0.124 | 0.88 |
| F vt1e2 rb1e5 | SimFilt | 1e2 | 1e5 | 0.12 | 0.153 / 12.37 | 64 | 0.112 | 0.75 |
| F vt0 rb1e6 | FAST | 0 | 1e6 | 0.04 | 0.193 / 12.56 | 81 | 0.205 | 1.13 |
| F vt0 rb1e5 | FAST | 0 | 1e5 | 0.04 | 0.106 / 12.43 | 51 | 0.315 | 6.2 |
| F vt0 rb1e5 | FAST | 0 | 1e5 | 0.12 | 0.315 / 13.50 | 132 | 0.243 | 5.2 |
| F vt1e2 rb1e5 | FAST | 1e2 | 1e5 | 0.12 | 0.307 / 13.49 | 128 | 0.244 | 5.7 |

- On the Simulink plant the blade flap model allows the lower pitch weight
  (stable even with DTc = 0.04) and beats the PI in speed, power and tower
  (vt 1e2, DTc 0.12).
- On FAST the tower oscillation of batch U is gone (tower 0.24-0.31 instead
  of 0.41-0.61), but the low pitch weight now excites other modes: with
  DTc = 0.04 the pitch rate variance is at the drivetrain mode (1.4-2.4 Hz)
  and 3P (0.4-0.8 Hz); with DTc = 0.12 at 2.4-6 Hz (step-wise pitch
  reference every 0.12 s).
- Next: notch filters for the drivetrain mode and 3P in the measured rotor
  speed, and a smooth (interpolated) pitch reference for DTc = 0.12.

## Batch N: notch filters and smooth pitch reference (flap model, FAST 18 m/s)

Flap model, omega weights, vt = 0, rMPC = [1 1e5], pMPC = 10; notch at
[0.605 1.64] Hz (3P at rated speed, FAST drivetrain mode), zetaNotchCtrl =
0.1 (depth 0.1, phase at the tower mode -10 deg).

| Name | DTc | Notch | Smooth | RotSpd std / max | Power std | Tower FA | Pitch rate std |
|---|---|---|---|---|---|---|---|
| PI | | | | 0.286 / 12.89 | 120 | 0.190 | 0.62 |
| U_FAST_omega_noTower (candidate omega) | 0.04 | | | 0.144 / 12.44 | 61 | 0.216 | 1.24 |
| F_FAST_vt0_rb1e5_DTc04 (candidate flap) | 0.04 | | | 0.106 / 12.43 | 51 | 0.315 | 6.2 |
| N_FAST_rb1e5_DTc04_notch | 0.04 | x | | 0.086 / 12.29 | 38 | 0.290 | 3.6 |
| N_FAST_rb1e5_DTc04_notch_smooth | 0.04 | x | x | 0.087 / 12.30 | 39 | 0.298 | 3.2 |
| N_FAST_rb1e5_DTc12_smooth | 0.12 | | x | 0.150 / 12.42 | 63 | 0.244 | 1.8 |
| N_FAST_rb1e5_DTc12_notch_smooth | 0.12 | x | x | 0.155 / 12.40 | 66 | 0.233 | 1.65 |

Variance per band (pitch rate [deg/s]^2 and tower acc. [m/s^2]^2 in
0-0.25 / 0.25-0.4 / 0.4-0.8 / 0.8-1.4 / 1.4-2.4 / 2.4-6.25 Hz):

| Name | Pitch rate | Tower acc. |
|---|---|---|
| PI | | 0.0010 / 0.0074 / 0.0184 / 0.0015 / 0.0007 / 0.0073 |
| omega candidate | 0.43 / 0.02 / 0.27 / 0.01 / 0.15 / 0.53 | 0.0016 / 0.0111 / 0.0248 / 0.0015 / 0.0007 / 0.0069 |
| flap candidate | 0.51 / 0.25 / 8.35 / 0.75 / 18.2 / 6.5 | 0.0026 / 0.0296 / 0.0547 / 0.0023 / 0.0016 / 0.0067 |
| DTc04 notch | 0.54 / 0.33 / 2.01 / 0.45 / 0.61 / 7.6 | 0.0036 / 0.0335 / 0.0356 / 0.0022 / 0.0007 / 0.0065 |
| DTc12 notch smooth | 0.47 / 0.09 / 0.31 / 0.06 / 0.11 / 1.5 | 0.0020 / 0.0163 / 0.0241 / 0.0015 / 0.0007 / 0.0071 |

- The notch filters remove the drivetrain-mode pitch activity (18.2 -> 0.6)
  and most of the 3P activity (8.4 -> 2.0). With DTc = 0.04 the rest is at
  2.4-6 Hz and the tower mode (0.25-0.4 Hz) stays 4.5 times the PI value.
- The smooth reference removes most of the 2.4-6 Hz content with DTc = 0.12
  (pitch rate 5.2 -> 1.8 deg/s).
- Best of this batch: DTc 0.12, notch and smooth reference (speed and power
  scatter about half of the PI, pitch rate 2.7 times, tower FA +23 %; the
  tower-mode band is still 2.2 times the PI value).
