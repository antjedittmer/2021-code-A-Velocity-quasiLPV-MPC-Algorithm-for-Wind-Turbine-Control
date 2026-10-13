function kf = flap_kalman(a, rhoP, DT, lambdaVec, pitchVec, cQ, cT, wecs, sigmaV, sigmaA)
% flap_kalman estimates the tower fore-aft and blade flap states of the
% qLPV-MPC model with blade flap (ABqL_WT_flap) with a Kalman filter.
%
%   z = [y_t, delta, v_t, ddelta] (deviations from the slow operating point)
%   dz/dt = Az*z + Bz*[dOmegaR; dbeta; dV] + Bz(:,3)*w   (w: wind error)
%   a = Az(3,:)*z + Bz(3,:)*[dOmegaR; dbeta; dV] + e    (measured tower acc.)
%
% Az, Bz are the tower/blade rows and columns of ABqL_WT_flap at the
% scheduling point; dOmegaR, dbeta and dV are the deviations of the measured
% rotor speed, pitch and wind from their 0.02 Hz means. The process noise is
% the difference between the measured hub wind and the rotor-effective wind
% (standard deviation sigmaV [m/s] per sample), the measurement noise sigmaA
% [m/s^2].
%
% Output: kf = [dy_t; ddelta; dv_t; dddelta; v_t] (increments over one
% sample and the tower velocity)

persistent z P rhoMean uPrev
if isempty(z)
    z = zeros(4, 1);
    P = zeros(4);
    rhoMean = rhoP(:);
    uPrev = zeros(3, 1);
end

fcMean = 0.02; % Hz, slow mean of the scheduling parameters (as fcTowerCtrl)
lam = exp(-2*pi*fcMean*DT);
rhoMean = lam*rhoMean + (1 - lam)*rhoP(:);
u = rhoP(:) - rhoMean;

[A, B] = ABqL_WT_flap(rhoP, lambdaVec, pitchVec, cQ, cT, wecs);
iz = [6 7 8 9];                   % y_t, delta, v_t, ddelta
Az = A(iz, iz);
Bz = [A(iz, 2), A(iz, 5), B(iz, 1)]; % rotor speed, pitch, wind
C = Az(3, :);
D = Bz(3, :);

% Exact zero-order-hold discretization
Md = expm([Az, Bz; zeros(3, 7)]*DT);
Ad = Md(1:4, 1:4);
Bd = Md(1:4, 5:7);

% Prediction with the input of the previous sample
zPred = Ad*z + Bd*uPrev;
Qd = (Bd(:, 3)*Bd(:, 3)')*sigmaV^2;
PPred = Ad*P*Ad' + Qd;

% Update with the measured acceleration
S = C*PPred*C' + sigmaA^2;
K = PPred*C'/S;
zNew = zPred + K*(a - C*zPred - D*u);
P = (eye(4) - K*C)*PPred;
P = (P + P')/2;

kf = [zNew - z; zNew(3)];
z = zNew;
uPrev = u;
end
