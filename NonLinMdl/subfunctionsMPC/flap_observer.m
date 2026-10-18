function flapMeas = flap_observer(a, v, rhoP, DT, lambdaVec, pitchVec, cQ, cT, wecs)
% flap_observer estimates the blade flap state of the qLPV-MPC model with
% blade flap (ABqL_WT_flap), which is not measured. The blade equation of
% Mdl2 is driven by the measured tower fore-aft acceleration a, the tower
% velocity v (from the acceleration, tower_states) and the deviations of
% rotor speed, pitch and wind from their slow mean values:
%
%   N*mb*rb*a + N*mb*rb^2*dddelta = -N*Kb*rb^2*delta - N*Bb*rb^2*ddelta + rb*Ft
%   Ft = dFt/dOmegaR*dOmegaR + dFt/dbeta*dbeta + dFt/dVe*(kWindEff*dV - v - rb*ddelta)
%
% The aerodynamic damping (dFt/dVe*rb^2) makes this open-loop estimate
% stable. Exact zero-order-hold discretization at the sample time DT.
%
% Inputs: * a: tower fore-aft acceleration [m/s^2]
%         * v: tower fore-aft velocity [m/s]
%         * rhoP: scheduling parameters [OmegaR [rad/s], beta [rad], V [m/s]]
%         * DT: sample time [s]
%         * lambdaVec, pitchVec, cQ, cT: tables, wecs: turbine parameters
%
% Output: flapMeas = [ddelta_k - ddelta_{k-1} increments: d(delta); d(ddelta)]

persistent x rhoMean
if isempty(x)
    x = [0; 0];
    rhoMean = rhoP(:);
end

fcMean = 0.02; % Hz, slow mean of the scheduling parameters (as fcTowerCtrl)
lam = exp(-2*pi*fcMean*DT);
rhoMean = lam*rhoMean + (1 - lam)*rhoP(:);
dRho = rhoP(:) - rhoMean;

[~, ~, sens] = ABqL_WT_flap(rhoP, lambdaVec, pitchVec, cQ, cT, wecs);
kWind = 1;
if isfield(wecs, 'kWindEff'), kWind = wecs.kWindEff; end

rb = wecs.rb;
mbN = wecs.N*wecs.mb;
mB = mbN*rb^2;                               % blade modal inertia
kB = wecs.N*wecs.Kb*rb^2;                    % blade stiffness
bB = wecs.N*wecs.Bb*rb^2 + rb^2*sens(6);     % structural + aerodynamic damping
force = rb*(sens(4)*dRho(1) + sens(5)*dRho(2) + sens(6)*(kWind*dRho(3) - v)) - mbN*rb*a;

Ac = [0, 1; -kB/mB, -bB/mB];
Bc = [0; 1/mB];
Md = expm([Ac, Bc; 0, 0, 0]*DT);
xNew = Md(1:2, 1:2)*x + Md(1:2, 3)*force;

flapMeas = xNew - x;
x = xNew;
end
