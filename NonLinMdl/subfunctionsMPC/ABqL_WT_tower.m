function [Aeval,Beval] = ABqL_WT_tower(rhoP,Rotor_Lambda,Rotor_Pitch,Rotor_cQ,Rotor_cT,wecs)
% ABqL_WT_tower calculates the continuous system (A) and input (B) matrix
% of the quasiLPV model with drivetrain, actuators and tower fore-aft
% motion. It is Mdl2 (Mdl_BianchiOL, initModel5MWNREL) without the blade
% flap and the tower side-side states, linearized at rhoP:
%
%   x = [ThetaS, OmegaR, OmegaG, Tg, beta, x_t, v_t]
%   u = [V, Tg_ref, beta_ref]
%
%   Jr*dOmegaR  = -Ks*ThetaS - Bs*(OmegaR - OmegaG) + Tr(kWindEff*V - v_t, OmegaR, beta)
%   mtb*dv_t    = -Kt*x_t - Bt*v_t + Ft(kWindEff*V - v_t, OmegaR, beta)
%   dx_t        = v_t
%
% Tr = 0.5*rho*pi*Rr^3*cQ(lambda,beta)*Ve^2, Ft = 0.5*rho*pi*Rr^2*cT(lambda,beta)*Ve^2,
% lambda = OmegaR*Rr/Ve.
%
% Inputs: * rhoP: scheduling parameters [OmegaR [rad/s], beta [rad], V [m/s]]
%         * Rotor_Lambda, Rotor_Pitch: grid of the cQ and cT tables
%         * Rotor_cQ, Rotor_cT: torque and thrust coefficient tables
%         * wecs: turbine parameters

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
Ng = wecs.Ng;
Jr = wecs.Jr;
Jgr = wecs.Jg * Ng^2;   % generator inertia in the LSS
Ks = wecs.Ks;
Bs = wecs.Bs;
mtb = wecs.mtb;         % tower fore-aft modal mass including blades (Mdl2: M(1,1))
Kt = wecs.Kt;
Bt = wecs.Bt;
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
dTrdVe = 2*kQ*Vbar*Cqb + kQ*Vbar^2*dCqdlambdab*dlambdadV; % w.r.t. Ve (tower velocity: dVe/dv_t = -1)
dFtdomega = kT*dCtdlambdab*dlambdadomega*Vbar^2;
dFtdbeta = kT*dCtdbetab*Vbar^2;
dFtdVe = 2*kT*Vbar*Ctb + kT*Vbar^2*dCtdlambdab*dlambdadV;

%% System matrix: x = [ThetaS, OmegaR, OmegaG, Tg, beta, x_t, v_t]
Aeval = zeros(7);
Aeval(1,2) = 1;
Aeval(1,3) = -1;
Aeval(2,:) = [-Ks, -Bs + dTrdomega, Bs, 0, dTrdbeta, 0, -dTrdVe]/Jr;
Aeval(3,:) = [Ks, Bs, -Bs, -Ng, 0, 0, 0]/Jgr;
Aeval(4,4) = -1/kappa_g;
Aeval(5,5) = -1/tau_g;
Aeval(6,7) = 1;
Aeval(7,:) = [0, dFtdomega, 0, 0, dFtdbeta, -Kt, -Bt - dFtdVe]/mtb;

%% Input matrix: u = [V, Tg_ref, beta_ref]
Beval = zeros(7,3);
Beval(2,1) = kWind*dTrdVe/Jr;   % hub height wind input: dVe/dV = kWind
Beval(4,2) = 1/kappa_g;
Beval(5,3) = 1/tau_g;
Beval(7,1) = kWind*dFtdVe/mtb;

end
