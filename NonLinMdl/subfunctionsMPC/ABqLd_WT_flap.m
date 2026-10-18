function [Av,Bv] = ABqLd_WT_flap(DT,rhoP,Rotor_Lambda,Rotor_Pitch,Rotor_cQ,Rotor_cT,wecs)
% ABqLd_WT_flap calculates the augmented discrete velocity-form model of
% the quasiLPV model with tower fore-aft motion and blade flap
% (ABqL_WT_flap).
%
% Inputs: * DT: sampling time
%         * rhoP: scheduling parameters [OmegaR, beta, V]
%
% Outputs: * Av, Bv: augmented system and input matrix
%          State order: [Tg, beta, OmegaR, v_t, dThetaS, dOmegaR, dOmegaG,
%          dTg, dBeta, dy_t, ddelta, dv_t, dddelta, dV, dTg_ref(k-1),
%          dBeta_ref(k-1)]

[Arho,Brho] = ABqL_WT_flap(rhoP,Rotor_Lambda,Rotor_Pitch,Rotor_cQ,Rotor_cT,wecs);

%% Discrete model (zero-order hold, exact)
nxc = size(Arho,1);
nuc = size(Brho,2);
Md = expm([Arho, Brho; zeros(nuc, nxc + nuc)]*DT);
Ak = Md(1:nxc, 1:nxc);
Bk = Md(1:nxc, nxc+1:end);

%% Velocity form augmented with the outputs y = [Tg, beta, OmegaR, v_t]
nyvec = [4, 5, 2, 8]; % in x = [ThetaS, OmegaR, OmegaG, Tg, beta, y_t, delta, v_t, ddelta]
ny = length(nyvec);
nx = nxc;
nd = 1;               % B column of the wind disturbance
nuvec = 2:3;          % B columns of the control inputs
nu = length(nuvec);

A = [eye(ny), Ak(nyvec,:);
    zeros(nx,ny), Ak];
B = [Bk(nyvec,nuvec); Bk(:,nuvec)];
Bv1 = [Bk(nyvec,nd); Bk(:,nd)];

% Augmented state [y; dx; dV; dU_{k-1}]: the actuators receive the
% reference with a one-sample delay, so the control increment acts on the
% plant one step later.
Av = [A, Bv1, B;
    zeros(nd, nx+ny), 0*eye(nd), zeros(nd, nu);
    zeros(nu, nx+ny+nd+nu)];
Bv = [zeros(nx+ny+nd, nu); eye(nu)];

end
