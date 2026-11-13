clc; clear; close all

%% Physical System (2nd order)
DT = 10^-2;
k = 0.5; b = 1; m = 2; f = 3;
Arho = [0, 1; -k/m, -b/m];
Brho = [0; f/m];

%% Discrete Plant for Delta x
Apre = (1/24)*DT^4*Arho^3 + (1/6)*DT^3*Arho^2 + 0.5*DT^2*Arho + DT*eye(2);
Ak   = Apre*Arho + eye(2);
Bk   = Apre*Brho;

%% Velocity Form
Av = [1,  Ak(1,:);
      zeros(2,1), Ak];
Bv = [Bk(1,1); Bk];

nx = 3; ni = 1; 
N = 18;
%N         = 10;     % Longer horizon = hard constraint has enough lookahead to brake
%N         = 40;     % Longer horizon for braking lookahead  


%% Prediction Matrices
L = zeros(nx*N, nx);
S = zeros(nx*N, ni*N);
L(1:nx, 1:nx) = Av;
S(1:nx, 1:ni) = Bv;
for idx = 1:N-1
    curr = idx*nx + (1:nx);
    prev = (idx-1)*nx + (1:nx);
    L(curr, :)            = Av * L(prev, :);
    S(curr, 1:(idx+1)*ni) = [Av * S(prev, 1:idx*ni), Bv];
end

%% Shared Settings
steps     = 5500;
pos_limit = 1.70;   % Hard ceiling
neg_limit = -0.5;
dU_max    = 0.011;
dU_min    = -dU_max;
ref       = 1.5;
R_val     = 350;

pos_limit = 1.55;   % Tight — only 0.05 above ref, forces early braking
R_val     = 300;     % Low R = aggressive tracking = large overshoot without constraints

R_val     = 500;    % Back to original — gives stable response

pos_limit = 1.55;   % Tight limit, only 0.05 above ref
ref       = 1.5;
dU_max    = 0.01;



% epsilon_buffer: tightens the constraint INWARD so the soft limit
% is slightly inside the hard limit, giving slack room to breathe.
% It must NOT be added to pos_limit — that would relax it outward.
epsilon_buffer = 0.03;   % Soft limit = pos_limit - epsilon_buffer = 1.52

Q = zeros(nx*N);
Q(1:(N-1)*nx, 1:(N-1)*nx)           = kron(eye(N-1), diag([1, 0, 0]));
Q((N-1)*nx+1:N*nx, (N-1)*nx+1:N*nx) = diag([1, 0, 0]);
R = kron(eye(N), R_val);

H = 2*(S'*Q*S + R);
H = (H + H') / 2;

T    = kron(eye(N), [1, 0, 0]);
opts = optimoptions('quadprog', 'Display', 'off');

Ac_u = [eye(N); -eye(N)];
bc_u = [dU_max*ones(N,1); -dU_min*ones(N,1)];

%% -----------------------------------------------------------------------
%% RUN 1: No state constraints
%% -----------------------------------------------------------------------
zk            = [0; 0; 0];
dU_guess      = zeros(N, 1);
hist_pos_free = zeros(steps, 1);
hist_du_free  = zeros(steps, 1);

for step = 1:steps
    g = 2 * S' * Q * (L * zk - ref);
    [dU_vec, ~, ~] = quadprog(H, g, Ac_u, bc_u, [], [], [], [], dU_guess, opts);

    if ~isempty(dU_vec)
        dU       = max(dU_min, min(dU_max, dU_vec(1)));
        dU_guess = [dU_vec(2:end); 0];
    else
        dU       = max(dU_min, min(dU_max, -zk(3)*0.1));
        dU_guess = zeros(N, 1);
    end

    zk = Av*zk + Bv*dU;
    hist_pos_free(step) = zk(1);
    hist_du_free(step)  = dU;
end

%% -----------------------------------------------------------------------
%% RUN 2: Soft state constraints (corrected epsilon direction)
%%
%% The soft constraint is:  T*S*dU - epsilon <= (pos_limit - epsilon_buffer) - T*L*zk
%% i.e. predicted positions must stay below (pos_limit - epsilon_buffer),
%% but can exceed it by at most epsilon (the slack variable, bounded above
%% by epsilon_buffer). This means the HARD ceiling is still pos_limit.

%% -----------------------------------------------------------------------
rho_slack = R_val * 1e4;
H_slack   = blkdiag(H, rho_slack);
H_slack   = (H_slack + H_slack') / 2;

zk            = [0; 0; 0];
dU_guess      = zeros(N+1, 1);
hist_pos_soft = zeros(steps, 1);
hist_du_soft  = zeros(steps, 1);

for step = 1:steps
    g       = 2 * S' * Q * (L * zk - ref);
    g_slack = [g; 0];

    TLzk = T * L * zk;

    % BUG FIX: subtract epsilon_buffer to tighten inward, not add it
    soft_pos_limit = pos_limit - epsilon_buffer;   % = 1.52
    soft_neg_limit = neg_limit + epsilon_buffer;   % = -0.47

    bc_s = [max(0, soft_pos_limit*ones(N,1) - TLzk);
            max(0, -soft_neg_limit*ones(N,1) + TLzk)];

    Ac_total = [Ac_u,         zeros(2*N, 1);
                [T*S; -T*S], -ones(2*N,  1)];
    bc_total = [bc_u; bc_s];

    % Slack bounded above by epsilon_buffer: beyond that, the hard limit
    % pos_limit would be breached, so we do not permit more slack
    lb = [-inf(N,1); 0             ];
    ub = [ inf(N,1); epsilon_buffer];

    [dU_vec, ~, ~] = quadprog(H_slack, g_slack, Ac_total, bc_total, ...
                               [], [], lb, ub, dU_guess, opts);

    if ~isempty(dU_vec)
        dU       = max(dU_min, min(dU_max, dU_vec(1)));
        dU_guess = [dU_vec(2:N); 0; dU_vec(N+1)];
    else
        dU       = max(dU_min, min(dU_max, -zk(3)*0.1));
        dU_guess = zeros(N+1, 1);
    end

    zk = Av*zk + Bv*dU;
    hist_pos_soft(step) = zk(1);
    hist_du_soft(step)  = dU;
end

%% -----------------------------------------------------------------------
%% RUN 3: Hard state constraints (no slack)
%% -----------------------------------------------------------------------
zk            = [0; 0; 0];
dU_guess      = zeros(N, 1);
hist_pos_hard = zeros(steps, 1);
hist_du_hard  = zeros(steps, 1);
hist_feas     = zeros(steps, 1);

for step = 1:steps
    g    = 2 * S' * Q * (L * zk - ref);
    TLzk = T * L * zk;

    Ac_pos  = T*S;
    bc_pos  = pos_limit*ones(N,1) - TLzk;
    Ac_neg  = -T*S;
    bc_neg  = neg_limit*(-ones(N,1)) + TLzk;

    Ac_full = [Ac_u; Ac_pos; Ac_neg];
    bc_full = [bc_u; bc_pos; bc_neg];

    [dU_vec, ~, exitflag] = quadprog(H, g, Ac_full, bc_full, ...
                                     [], [], [], [], dU_guess, opts);
    if exitflag == 1
        hist_feas(step) = 1;
        dU       = max(dU_min, min(dU_max, dU_vec(1)));
        dU_guess = [dU_vec(2:end); 0];
    else
        % Fallback: upper limit only
        [dU_vec2, ~, ef2] = quadprog(H, g, [Ac_u; Ac_pos], [bc_u; bc_pos], ...
                                     [], [], [], [], dU_guess, opts);
        if ef2 == 1
            hist_feas(step) = 2;
            dU       = max(dU_min, min(dU_max, dU_vec2(1)));
            dU_guess = [dU_vec2(2:end); 0];
        else
            % Fallback: rate limits only
            [dU_vec3, ~, ef3] = quadprog(H, g, Ac_u, bc_u, ...
                                         [], [], [], [], dU_guess, opts);
            if ef3 == 1
                hist_feas(step) = 3;
                dU       = max(dU_min, min(dU_max, dU_vec3(1)));
                dU_guess = [dU_vec3(2:end); 0];
            else
                hist_feas(step) = 4;
                dU       = max(dU_min, min(dU_max, -zk(3)*0.1));
                dU_guess = zeros(N, 1);
            end
        end
    end

    zk = Av*zk + Bv*dU;
    hist_pos_hard(step) = zk(1);
    hist_du_hard(step)  = dU;
end

fprintf('Hard constraint feasibility:\n');
fprintf('  Full (both limits):  %d steps\n', sum(hist_feas==1));
fprintf('  Upper limit only:    %d steps\n', sum(hist_feas==2));
fprintf('  Rate limits only:    %d steps\n', sum(hist_feas==3));
fprintf('  Emergency:           %d steps\n', sum(hist_feas==4));

%% -----------------------------------------------------------------------
%% Plot
%% -----------------------------------------------------------------------
time = (0:steps-1)*DT;

figure('Position', [100 100 960 750]);

subplot(3,1,1); hold on;
plot(time, hist_pos_free, 'b',   'LineWidth', 2);
plot(time, hist_pos_soft, 'r--', 'LineWidth', 2);
plot(time, hist_pos_hard, 'm:',  'LineWidth', 2.5);
yline(pos_limit,              'k-',  'LineWidth', 1.5);
yline(pos_limit-epsilon_buffer,'k:', 'LineWidth', 1.0);  % Soft limit line
yline(neg_limit,              'k-',  'LineWidth', 1.5);
yline(ref,                    'g--', 'LineWidth', 1.5);
legend('No constraints','Soft constraints','Hard constraints', ...
       'Hard limit','Soft limit','','Reference','Location','southeast');
title(sprintf('Position Comparison  (R = %g,  pos\\_limit = %.2f,  soft\\_limit = %.2f)', ...
      R_val, pos_limit, pos_limit-epsilon_buffer));
xlabel('Time (s)'); ylabel('Position'); grid on;

subplot(3,1,2); hold on;
stairs(time, hist_du_free, 'b',   'LineWidth', 2);
stairs(time, hist_du_soft, 'r--', 'LineWidth', 2);
stairs(time, hist_du_hard, 'm:',  'LineWidth', 2.5);
legend('No constraints','Soft constraints','Hard constraints','Location','northeast');
title('\Delta u(k) Comparison');
xlabel('Time (s)'); ylabel('\Delta u'); grid on;

subplot(3,1,3);
area(time, (hist_feas==1)*1,   'FaceColor',[0.2 0.7 0.2],'EdgeColor','none'); hold on;
area(time, (hist_feas==2)*0.7, 'FaceColor',[1.0 0.8 0.0],'EdgeColor','none');
area(time, (hist_feas==3)*0.4, 'FaceColor',[1.0 0.4 0.0],'EdgeColor','none');
area(time, (hist_feas==4)*0.1, 'FaceColor',[0.8 0.0 0.0],'EdgeColor','none');
legend('Full (both limits)','Upper only','Rate only','Emergency','Location','northeast');
title('Hard Constraint Feasibility');
xlabel('Time (s)'); ylim([0 1.2]); yticks([]); grid on;