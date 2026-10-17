function [beta, ipc] = ipc_dlr(betaColl, MyKNm, azDeg, genSpd, twrAcc, ipcPar, DT, maxRateDeg)
% ipc_dlr individual pitch control (IPC) on top of a collective pitch
% reference: DLR IPC PI controller with tower damping (IPC-Regler1 in
% Sandbox/filesCompareDifferentConstraints/WECSCtrl/MdLCtrl_IPC.slx, gains
% in subfunctionsUtil/initDLRparameter.m, thesis section 4.2), discretized
% with the Tustin method at the sample time DT.
%
%   f_i     = BP(IPC_K_1P*My_i), BP = c1*s/(d1*s^2 + d2*s + d3) (1P band-pass)
%   [Mc;Ms] = 2/3*[cos(psi_i); sin(psi_i)]*f      (multiblade coordinates)
%   q_c,s   = kScale*(KI/s + KP)*M_c,s             (PI, same gains)
%   q_i     = cos(psi_i + off)*q_c + sin(psi_i + off)*q_s
%   ipc_i   = g(omega_g)*(F(q_i) + towerOn*kp_Turm*LP(a_t)), F = 100/(s^2 + 20*s + 100)
%   beta_i  = betaColl + ipc_i, rate limit maxRateDeg, limits -2 to 90 deg
%
% psi_i = azimuth of blade i (psi_1 = azDeg in rad, + 120 deg per blade);
% g fades the IPC in between 80 and 90 % of the rated generator speed.
%
% Sign: with the FAST blade root moments RootMOoP and azimuth the IPC
% offsets are added to the collective pitch (in the DLR model they are
% subtracted); with the DLR sign the IPC is unstable on FAST (gain set 2)
% or increases the moments (sets 1 and 3). On FAST (18 m/s, flapNotch,
% 150 s): gain set 2 is too aggressive (blade root moment std +18 %),
% kScale = 0.3 reduces it by 22 %; sets 1 and 3 act hardly at all (-2 %);
% the tower term increases the tower acceleration (the qLMPC damps the
% tower), see experiments/EXPERIMENTS.md.
%
% Inputs: betaColl collective pitch [rad], MyKNm blade root out-of-plane
% moments [kNm] (3x1), azDeg rotor azimuth [deg], genSpd generator speed
% [rad/s], twrAcc tower fore-aft acceleration [m/s^2], ipcPar = [gainSet,
% kScale, towerOn] (gainSet 0: IPC off, 1 to 3: tuned gain sets of
% initDLRparameter switchIPC; kScale: factor of the PI gains; towerOn: 1 with
% the tower damping term kp_Turm), DT [s], maxRateDeg pitch rate limit [deg/s].
% Outputs: blade pitch commands beta [rad] (3x1) and IPC pitch offsets ipc
% [rad] (3x1, before the rate limit)

persistent bpZ actZ intC intS lpZ wFilt betaPrev
if isempty(bpZ)
    bpZ = zeros(2, 3);
    actZ = zeros(2, 3);
    intC = 0;
    intS = 0;
    lpZ = 0;
    wFilt = genSpd;
    betaPrev = betaColl*ones(3, 1);
end

ipc = zeros(3, 1);
gainSet = ipcPar(1);
kScale = ipcPar(2);
towerOn = ipcPar(3);
if gainSet < 1
    beta = betaColl*ones(3, 1);
    betaPrev = beta;
    return
end

% Parameters (initDLRparameter)
K1P = 1e-6;                  % moment scaling (Nm -> MNm)
c1 = 2*pi/4.95; d1 = 1; d3 = (2*pi/4.95)^2;
offset1P = 25*pi/180;        % azimuth offset of the back-transformation
TTwr = 0.1304*pi;            % time constant of the tower acceleration low-pass
vGenRated = 122.9096;        % rated generator speed [rad/s]
switch gainSet
    case 1
        KI = 1.0495e-05; KP = 0.00042667; d2 = 1.4983; kpTwr = 1e-05;
    case 2
        KI = 0.00028; KP = 0.02786; d2 = 0.3120; kpTwr = 0.0086982;
    otherwise
        KI = 4.45e-05; KP = -0.00010109; d2 = 1.925; kpTwr = 1.9362e-05;
end

% Fade-in with the filtered generator speed (first-order low-pass 0.25 Hz)
alpha = 1 - exp(-2*pi*0.25*DT);
wFilt = wFilt + alpha*(genSpd - wFilt);
g = min(max((wFilt/vGenRated - 0.8)/0.1, 0), 1);

% Band-pass of the scaled blade root moments
[bBP, aBP] = tustin2([0 c1 0], [d1 d2 d3], DT);
f = zeros(3, 1);
for k = 1:3
    [f(k), bpZ(:, k)] = biquad(bBP, aBP, K1P*1000*MyKNm(k), bpZ(:, k));
end

% Multiblade coordinate transformation and PI controllers
psi = azDeg*pi/180 + [0; 2*pi/3; 4*pi/3];
Mc = 2/3*(cos(psi)'*f);
Ms = 2/3*(sin(psi)'*f);
intC = intC + Mc*DT;
intS = intS + Ms*DT;
qc = kScale*(KI*intC + KP*Mc);
qs = kScale*(KI*intS + KP*Ms);
q = cos(psi + offset1P)*qc + sin(psi + offset1P)*qs;

% Tower damping (low-pass of the tower acceleration)
[bLP, aLP] = tustin1(1, [TTwr 1], DT);
yLP = bLP(1)*twrAcc + lpZ;
lpZ = bLP(2)*twrAcc - aLP(2)*yLP;
twr = towerOn*kpTwr*yLP;

% Actuator filter of the cyclic signals, fade-in, blade pitch commands
[bAct, aAct] = tustin2([0 0 100], [1 20 100], DT);
for k = 1:3
    [qF, actZ(:, k)] = biquad(bAct, aAct, q(k), actZ(:, k));
    ipc(k) = g*(qF + twr);
end
betaRaw = min(max(betaColl + ipc, -2*pi/180), 90*pi/180); % FAST sign convention
dMax = maxRateDeg*pi/180*DT;
beta = betaPrev + min(max(betaRaw - betaPrev, -dMax), dMax);
betaPrev = beta;
end

function [b, a] = tustin2(num, den, DT)
% Tustin discretization of (num(1)*s^2 + num(2)*s + num(3))/(den(1)*s^2 + den(2)*s + den(3))
K = 2/DT;
b = [num(1)*K^2 + num(2)*K + num(3), 2*(num(3) - num(1)*K^2), num(1)*K^2 - num(2)*K + num(3)];
a = [den(1)*K^2 + den(2)*K + den(3), 2*(den(3) - den(1)*K^2), den(1)*K^2 - den(2)*K + den(3)];
b = b/a(1);
a = a/a(1);
end

function [b, a] = tustin1(num, den, DT)
% Tustin discretization of num/(den(1)*s + den(2))
K = 2/DT;
b = [num, num];
a = [den(1)*K + den(2), den(2) - den(1)*K];
b = b/a(1);
a = a/a(1);
end

function [y, z] = biquad(b, a, x, z)
% Direct form II transposed
y = b(1)*x + z(1);
z = [b(2)*x - a(2)*y + z(2); b(3)*x - a(3)*y];
end
