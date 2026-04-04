clc; clear; close all

% Physical System (2nd order)
DT = 10^-2;
k = 0.5; b = 1; m = 2; f = 3;
Arho = [0, 1; -k/m, -b/m];
Brho = [0; f/m];

% Discrete Plant for Delta x
Apre = (1/24)*DT^4*Arho^3 + (1/6)*DT^3*Arho^2 + 0.5*DT^2*Arho + DT*eye(2);
Ak = Apre*Arho + eye(2);
Bk = Apre*Brho;

% YOUR VELOCITY FORM DEFINITION:
% z(k) = [x_absolute; delta_x1; delta_x2]
% x_abs(k+1) = x_abs(k) + delta_x1(k+1)
%            = x_abs(k) + Ak(1,:) * delta_x(k) + Bk(1) * delta_u(k)

ny = 1; % One absolute state to track/limit
Av = [1,  Ak(1,:);
    zeros(2,1), Ak];

Bv = [Bk(1,1);
    Bk];

nx = 3; ni = 1; N = 18;

% Re-implementing your exact recursive L and S logic
L = zeros(nx*N, nx);
S = zeros(nx*N, ni*N);
L(1:nx, 1:nx) = Av;
S(1:nx, 1:ni) = Bv;

for idx = 1:N-1
    curr = idx*nx + (1:nx);
    prev = (idx-1)*nx + (1:nx);
    L(curr, :) = Av * L(prev, :);
    S(curr, 1:(idx+1)*ni) = [Av * S(prev, 1:idx*ni), Bv];
end

%% Test
% t = 0: DT : 100;
% u = zeros(size(t)); u(1) = 1;
% sys = ss(Av,Bv, [1,0,0],0,DT);
% lsim(sys,u,t)

%% Simulation Setup
steps = 3500;
zk = [0; 0; 0]; % [pos=0, delta_pos=0, delta_vel=0]
pos_limit = 1.7; %1.7;
neg_limit = -0.5;
dU_max = 0.01;
dU_min = - dU_max;

% Q targets the first state (absolute position)
Q = zeros(nx*N);
Q(1: (N-1)*nx, 1: (N-1)*nx) = kron(eye(N-1), diag([1, 0, 0]));
Q((N-1)*nx + 1 :N*nx, (N-1)*nx + 1 :N*nx ) = diag([1, 0, 0]);

R = zeros(N);
R(1:N-1, 1:N-1) = kron(eye(N-1), 500);
R(N, N) = 500;
H = 2 * (S' * Q * S + R);
H = (H+H')/2;

% --- Setup before loop ---
rho_slack = 1e9; % Penalty for boundary violation
H_slack = blkdiag(H, rho_slack); % Expand H to handle [dU_vec; epsilon]
H_slack = (H_slack + H_slack') / 2;

ref = 1.5;
ref_vec =  ref * ones(nx*N, 1); % Target 5, should be stopped at 2

% Selector for the absolute position (state 1)
T = kron(eye(N), [1, 0, 0]);

history_pos = zeros(steps,1);
history_du = zeros(steps,1);

dU_guess = zeros(N, 1);

for k = 1:steps
    g = 2 * S' * Q * (L * zk - ref);
    g = 2 * S' * Q * (L * zk - ref);
    g_slack = [g; 0];

    % Constraints
    % 2. Hard Rate Limits (No slack column)
    Ac_u = [eye(N); -eye(N)];
    bc_u = [dU_max * ones(N,1); -dU_min * ones(N,1)];


    % 3. Soft State Limits (Adding slack variable column)
    % Logic: x <= limit + eps  =>  (T*S)*dU - eps <= limit - T*L*zk
    Ac_s = [T * S; -T * S];
    epsilon_buffer = 0.05;
    bc_s = [(pos_limit + epsilon_buffer) * ones(N, 1) - T * L * zk;
        (-neg_limit + epsilon_buffer) * ones(N, 1) + T * L * zk];

    
    % --- THE CRITICAL CHANGE ---
    % We add a column of -1 to state constraints and 0 to rate constraints
    Ac_total = [Ac_u,  zeros(2*N, 1);  % Rate limits don't bend
                Ac_s, -ones(2*N, 1)]; % State limits DO bend
    bc_total = [bc_u; bc_s];
    
    % Lower bound: Slack (the last variable) must be >= 0
    lb = [-inf(N, 1); 0]; 
    % Initialize dU_guess before the loop: dU_guess = zeros(N, 1);

    
    [dU_vec, ~, exitflag] = quadprog(H_slack, g_slack, Ac_total, bc_total, [], [],  lb, [], dU_guess, optimoptions('quadprog','Display','off'));
    %[dU_vec, ~, exitflag] = quadprog(H, g, Ac_u, bc_u, [], [], [], [], [], optimoptions('quadprog','Display','off'));

    if ~isempty(dU_vec)
        dU = max(dU_min, min(dU_max, dU_vec(1))); %dU_vec(1);
        dU_guess = [dU_vec(2:end); 0];
    else
        [dU_vec, ~, exitflag] = quadprog(H_slack, g_slack, Ac_u, bc_u, [], [],  lb, [], dU_guess, optimoptions('quadprog','Display','off'));

        if ~isempty(dU_vec)
        dU = max(dU_min, min(dU_max, dU_vec(1))); %dU_vec(1);
        dU_guess = [dU_vec(2:end); 0];
        else
        

        % Emergency braking if everything fails
        backup_val = -zk(2) * 0.001; % Proportional braking based on velocity
        dU = max(dU_min, min(dU_max, backup_val));
        end
    end

    zk = Av * zk + Bv * dU;
    history_pos(k) = zk(1);
    history_du(k) = dU;

    dU_vec_last{k} = dU_vec;
end

% Plotting
figure;
time = 0:DT: (steps-1)*DT;
subplot(2,1,1); plot(time, history_pos, 'LineWidth', 2); yline(pos_limit, 'r--','linewidth',2); yline(ref, 'g--','linewidth',2); title('Absolute Position x(k)'); grid on;
subplot(2,1,2); stairs(time,history_du, 'LineWidth', 2); title('\Delta u(k)'); grid on;