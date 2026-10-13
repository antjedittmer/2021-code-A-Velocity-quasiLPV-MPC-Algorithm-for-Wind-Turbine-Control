function [Aeval,Beval,sens] = ABqL_WT_flap(rhoP,Rotor_Lambda,Rotor_Pitch,Rotor_cQ,Rotor_cT,wecs)
% ABqL_WT_flap calculates the continuous system (A) and input (B) matrix
% of the quasiLPV model with drivetrain, actuators, tower fore-aft motion
% and blade flap. It is Mdl2 (Mdl_BianchiOL, initModel5MWNREL) without the
% tower side-side state, linearized at rhoP:
%
%   x = [ThetaS, OmegaR, OmegaG, Tg, beta, y_t, delta, v_t, ddelta]
%   u = [V, Tg_ref, beta_ref]
%
%   Jr*dOmegaR = -Ks*ThetaS - Bs*(OmegaR - OmegaG) + Tr(Ve, OmegaR, beta)
%   Mfa*[dv_t; dddelta] = -Kfa*[y_t; delta] - Bfa*[v_t; ddelta] + [1; rb]*Ft(Ve, OmegaR, beta)
%   Mfa = [mtb, N*mb*rb; N*mb*rb, N*mb*rb^2], Kfa = diag(Kt, N*Kb*rb^2),
%   Bfa = diag(Bt, N*Bb*rb^2)
%
% with the relative rotor-effective wind Ve = kWindEff*V - v_t - rb*ddelta,
% Tr = 0.5*rho*pi*Rr^3*cQ(lambda,beta)*Ve^2, Ft = 0.5*rho*pi*Rr^2*cT(lambda,beta)*Ve^2,
% lambda = OmegaR*Rr/Ve.
%
% Inputs: * rhoP: scheduling parameters [OmegaR [rad/s], beta [rad], V [m/s]]
%         * Rotor_Lambda, Rotor_Pitch: grid of the cQ and cT tables
%         * Rotor_cQ, Rotor_cT: torque and thrust coefficient tables
%         * wecs: turbine parameters
%
% Output sens: aerodynamic sensitivities [dTr/dOmegaR, dTr/dbeta, dTr/dVe,
% dFt/dOmegaR, dFt/dbeta, dFt/dVe] (for the blade flap observer)

persistent lambdaVec pitchVec cQ cT Cqdlambda Cqdbeta Ctdlambda Ctdbeta

if isempty(lambdaVec)
    lambdaVec = Rotor_Lambda;
    pitchVec = Rotor_Pitch;
    cQ = Rotor_cQ;
    cT = Rotor_cT;
    [Cqdlambda,Cqdbeta] = gradient(Rotor_cQ,mean(diff(Rotor_Lambda)),mean(diff(Rotor_Pitch)));
    [Ctdlambda,Ctdbeta] = gradient(Rotor_cT,mean(diff(Rotor_Lambda)),mean(diff(Rotor_Pitch)));
end

%% Operating point
omegabar = rhoP(1);
beta_bar = rhoP(2);
Vbar = rhoP(3);
% Rotor-effective wind for the aerodynamics, as in the WECS models:
% Ve = kWindEff*V (V: hub height wind)
kWind = 1;
if isfield(wecs, 'kWindEff'), kWind = wecs.kWindEff; end
Vbar = kWind*Vbar;

%% Physical and model constants
rho = 1.2250;
Rr = wecs.Rr;
rb = wecs.rb;           % aerodynamic center on the blade
Ng = wecs.Ng;
Jr = wecs.Jr;
Jgr = wecs.Jg * Ng^2;   % generator inertia in the LSS
Ks = wecs.Ks;
Bs = wecs.Bs;
mbN = wecs.N * wecs.mb; % modal mass of all blades
Mfa = [wecs.mtb, mbN*rb; mbN*rb, mbN*rb^2];
Kfa = diag([wecs.Kt, wecs.N*wecs.Kb*rb^2]);
Bfa = diag([wecs.Bt, wecs.N*wecs.Bb*rb^2]);
tau_g = 0.05;           % pitch actuator time constant, as tau in DD_Mdl_Actuator.sldd
kappa_g = 0.01;         % generator torque actuator time constant, as kappa in DD_Mdl_Actuator.sldd

%% Aerodynamic sensitivities
Lambda_bar = omegabar*Rr/Vbar;
dlambdadomega = Rr/Vbar;
dlambdadV = -omegabar*Rr/Vbar^2;
% stay within the table grid (interp2 returns NaN outside)
Lambda_bar = min(max(Lambda_bar, lambdaVec(1)), lambdaVec(end));
beta_bar = min(max(beta_bar, pitchVec(1)), pitchVec(end));

Cqb = interp2(lambdaVec,pitchVec,cQ,Lambda_bar,beta_bar);
dCqdlambdab = interp2(lambdaVec,pitchVec,Cqdlambda,Lambda_bar,beta_bar);
dCqdbetab = interp2(lambdaVec,pitchVec,Cqdbeta,Lambda_bar,beta_bar);
Ctb = interp2(lambdaVec,pitchVec,cT,Lambda_bar,beta_bar);
dCtdlambdab = interp2(lambdaVec,pitchVec,Ctdlambda,Lambda_bar,beta_bar);
dCtdbetab = interp2(lambdaVec,pitchVec,Ctdbeta,Lambda_bar,beta_bar);

kQ = 0.5*rho*pi*Rr^3;
kT = 0.5*rho*pi*Rr^2;
dTrdomega = kQ*dCqdlambdab*dlambdadomega*Vbar^2;
dTrdbeta = kQ*dCqdbetab*Vbar^2;
dTrdVe = 2*kQ*Vbar*Cqb + kQ*Vbar^2*dCqdlambdab*dlambdadV; % w.r.t. Ve (dVe/dv_t = -1, dVe/dddelta = -rb)
dFtdomega = kT*dCtdlambdab*dlambdadomega*Vbar^2;
dFtdbeta = kT*dCtdbetab*Vbar^2;
dFtdVe = 2*kT*Vbar*Ctb + kT*Vbar^2*dCtdlambdab*dlambdadV;

%% System matrix: x = [ThetaS, OmegaR, OmegaG, Tg, beta, y_t, delta, v_t, ddelta]
Aeval = zeros(9);
Aeval(1,2) = 1;
Aeval(1,3) = -1;
Aeval(2,:) = [-Ks, -Bs + dTrdomega, Bs, 0, dTrdbeta, 0, 0, -dTrdVe, -dTrdVe*rb]/Jr;
Aeval(3,:) = [Ks, Bs, -Bs, -Ng, 0, 0, 0, 0, 0]/Jgr;
Aeval(4,4) = -1/kappa_g;
Aeval(5,5) = -1/tau_g;
Aeval(6,8) = 1;
Aeval(7,9) = 1;
% generalized forces [1; rb]*Ft on the tower and blade coordinates
qF = [1; rb];
Aeval(8:9,:) = Mfa\([zeros(2,1), qF*dFtdomega, zeros(2,2), qF*dFtdbeta, -Kfa, ...
    -Bfa - qF*dFtdVe*[1, rb]]);

%% Input matrix: u = [V, Tg_ref, beta_ref]
Beval = zeros(9,3);
Beval(2,1) = kWind*dTrdVe/Jr;   % hub height wind input: dVe/dV = kWind
Beval(4,2) = 1/kappa_g;
Beval(5,3) = 1/tau_g;
Beval(8:9,1) = Mfa\(qF*kWind*dFtdVe);

sens = [dTrdomega, dTrdbeta, dTrdVe, dFtdomega, dFtdbeta, dFtdVe];

end
