function [Av,Bv] = ABqLd_WT_beta(DT,rhoP,Rotor_Lambda,Rotor_Pitch,Rotor_cQ,wecs)
% ABqLd calculates system (A) and input (B) matrix for quasiLPV
% representation from ABgL and augments the states with integrating
% disturbances.
%
% Inputs: * DT: Sampling time
%         * rho: Scheduling parameter vector: rho = V
%
% Outputs: * A: Augmented system matrix A(rho) at current rho
%          * B: Augmented Input matrix B(rho) at current rho
%          State order: [Tg, beta, OmegaR, dThetaS, dOmegaR, dOmegaG, dTg,
%          dBeta, dV, dTg_ref(k-1), dBeta_ref(k-1)]

% Pablo S.G. Cisneros, Herbert Werner, ICS TUHH
% modified for WECS simulation: Antje Dittmer

% Calculate system and input matrices at current wind speed
[Arho,Brho] = ABqL_WT_beta(rhoP,Rotor_Lambda,Rotor_Pitch,Rotor_cQ,wecs);

%% Calculate discrete forms of A,B (zero-order hold, exact)
% The 4th order Taylor series used before diverges for DT larger than the
% fast actuator time constants (kappa_g = 0.01 s, tau_g = 0.05 s), e.g.
% for a controller sample time DTc = 0.04 s.
nxc = size(Arho,1);
nuc = size(Brho,2);
Md = expm([Arho, Brho; zeros(nuc, nxc + nuc)]*DT);
Ak = Md(1:nxc, 1:nxc);
Bk = Md(1:nxc, nxc+1:end);

%% Calculate new matrices augmented system
ny = 3; % outputs y = [Tg, beta, OmegaR]: the absolute rotor speed is tracked
nx = size(Arho,1);
nyvec = [nx-1, nx, 2]; % Tg, beta, OmegaR in x = [ThetaS, OmegaR, OmegaG, Tg, beta]
nd = 1;% B: columns disturbance input
nuvec = 2:3; % B: columns control input
Cdiag = eye(ny); %diag([97,1]); %eye(ny); %diag([97,1]);


nu = length(nuvec);

A = [eye(ny) Cdiag * Ak(nyvec,:);
    zeros(size(Ak,1),ny), Ak];
B = [Cdiag * Bk(nyvec,nuvec); Bk(:,nuvec)];
Bv1 = [Cdiag * Bk(nyvec,nd); Bk(:,nd)];

% Augmented state [y; dx; dV; dU_{k-1}]: the actuators receive the
% reference with a one-sample delay (Delay blocks in Mdl_Actuator), so the
% control increment acts on the plant one step later.
Av = [A, Bv1, B;
    zeros(nd, nx+ny), 0*eye(nd), zeros(nd, nu);
    zeros(nu, nx+ny+nd+nu)];
Bv = [zeros(nx+ny+nd, nu); eye(nu)];


end