clear; clc; close all;
% addpath(genpath('C:\Users\CompuTop\seadrive_root\Emmanuel\Meine Bibliotheken\My Library\MATLAB-Libs\sedumi-master'))
addpath(genpath('C:\Users\CompuTop\seadrive_root\Emmanuel\Meine Bibliotheken\My Library\MATLAB-Libs\YALMIP-master'))
addpath(genpath("C:\Users\CompuTop\seadrive_root\Emmanuel\Meine Bibliotheken\My Library\MATLAB-Libs\sedumi"))
addpath(genpath("C:\Users\CompuTop\seadrive_root\Emmanuel\Meine Bibliotheken\My Library\MATLAB-Libs\matlab2tikz-master"))
load('data_problem2.mat');

% ---------------------------------------------------------------
% Load the TRUE A and B matrices from the provided function
% ---------------------------------------------------------------
A_true = [1, -18.4542,-3.3219, 2.2083, 0.1722;
          0, 0.6807, 0.2618, 0, 0;
          0, 0, 0.6807, 0, 0;
          0, 0, 0, 0.8519, 0.1365;
          0, 0, 0, 0, 0.8519];

B_true = [-0.4539, 0.0094;
           0.0575, 0;
           0.3193, 0;
           0, 0.01155;
           0, 0.1481];

%% ---------------------------------------------------------------
% Loop through datasets 1, 2, and 3
% ---------------------------------------------------------------
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 2(a)-(b): Data informativity for controllability and stabilizability\n');
fprintf('--------------------------------------------------------------\n');
fprintf('\n');


for d = 1:3

    % select dataset
    x_off = eval(sprintf('x_off_%d', d));   % 5 x (N+1)
    u_off = eval(sprintf('u_off_%d', d));   % 2 x (N+1)

    [n, Np1] = size(x_off);   % n states, N+1 samples
    [m, ~   ] = size(u_off);  % m inputs
    N = Np1 - 1;

    % format data with time in rows for H(.)
    x_t = x_off.';            % (N+1) x 5
    u_t = u_off.';            % (N+1) x 2

    % build X-, U-, X+
    Hx  = H(1, x_t(1:N,   :));   % X_-
    Hu  = H(1, u_t(1:N,   :));   % U_-
    Hxp = H(1, x_t(2:N+1, :));   % X_+

    if d == 3
        Hx_3  = Hx;
        Hu_3  = Hu;
        Hxp_3 = Hxp;
    end

    % check informativity
    fprintf('\n--------------------------------------------------------------\n');
    fprintf(' Testing dataset %d for informativity for system identification\n', d);
    fprintf('--------------------------------------------------------------\n\n');

    Hxu = [Hx; Hu];   % (n+m) x N

    if rank(Hxu) == (n + m)
        fprintf('Data %d is informative for system identification.\n', d);

        % identification of [A  B]
        M = Hxp * pinv(Hxu);    % size = 5 x (5+2) = 5 x 7

        A_id = M(:, 1:n);
        B_id = M(:, n+1:end);

        fprintf('Identified A_%d = \n'); disp(A_id);
        fprintf('Identified B_%d = \n'); disp(B_id);

        % compute identification errors
        A_error = A_id - A_true;
        B_error = B_id - B_true;

        fprintf('A error for dataset %d = \n'); disp(A_error);
        fprintf('B error for dataset %d = \n'); disp(B_error);

    else
        fprintf('Data %d is NOT informative for system identification.\n', d);
    end

    fprintf('\n');
    fprintf('--------------------------------------------------------------\n');
    fprintf(' Testing data %d for informativitiy for controllability\n', d);
    fprintf('--------------------------------------------------------------\n');
    fprintf('\n');

    ok = '';   % Flag for result

    % First condition: rank(X_+) must be n
    if rank(Hxp) < n
        ok = ' NOT';
    else
        % Compute sigma = eigenvalues of (X_- * pinv(X_+))
        sigma = eig(Hx * pinv(Hxp));

        % Loop over lambdas
        for i = 1:n
            lambda = 1/sigma(i);
            fprintf('    Test %02d/%02d: lambda = %10.3e ... ', i, n, lambda);

            if rank(Hxp - lambda*Hx) == n
                fprintf('OK\n');
            else
                fprintf('NOT OK\n');
                ok = ' NOT';
                break
            end
        end
    end

    fprintf('\n');
    fprintf('Data is%s informative for controllability\n', ok);

    fprintf('\n');
    fprintf('--------------------------------------------------------------\n');
    fprintf(' Testing data %d for informativitiy for stabilizability\n', d);
    fprintf('--------------------------------------------------------------\n');
    fprintf('\n');

    ok = '';

    if rank(Hxp - Hx) < n
        ok = ' NOT';
    else
        sigma = eig(Hx * pinv(Hxp - Hx));

        for i = 1:n
            lambda = 1/sigma(i) + 1;
            fprintf('    Test %02d/%02d: lambda = %10.3e ... ', i, n, lambda);

            if rank(Hxp - lambda*Hx) == n
                fprintf('OK\n');
            else
                fprintf('NOT OK\n');
                ok = ' NOT';
                break
            end
        end
    end

    fprintf('\n');
    fprintf('Data is%s informative for stabilizability\n', ok);

end


%% ------------------------------------------------------------------
% Data informativity for stabilization by state feedback
% -------------------------------------------------------------------
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 2(c): Data informativity for stabilization by state feedback\n');
fprintf('--------------------------------------------------------------\n');
fprintf('\n');

% Solver Settings
sdpopt = sdpsettings('solver','sedumi','verbose',1);
sdpopt.sedumi.numtol = 1e-6;
sdpopt.sedumi.eps = 1e-9;

N = size(Hx_3,2);
n = size(Hx_3,1);

% Formulate the LMI problem
theta = sdpvar(N,n);            % Define theta as unknown variable in LMI
X     = sdpvar(n,n,'symmetric');    % Define X (placeholder for Hx*theta) as a symmetric matrix

F  = [ [ X Hxp_3*theta ; theta'*Hxp_3' X ] >= 1e-6 ];  % The main LMI
F  = [ F, X == Hx_3*theta ];                         % X = Hx*theta                              

st = optimize (F, [], sdpopt);

ths = value(theta);
Xs  = value(X);

% Stabilizing gain
K_stab = Hu_3 * ths / Xs;
fprintf('K_stab = \n');
disp(K_stab);

% Optional validation with true system matrices
eig_closed_loop = eig(A_true + B_true*K_stab);
fprintf('eig(A_true + B_true*K_stab) = \n');
disp(eig_closed_loop);


%% -------------------------------------------------------------
% Problem 2(d): Closed-loop system simulation
% --------------------------------------------------------------
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 2(d): Closed-loop system simulation\n');
fprintf('--------------------------------------------------------------\n');
fprintf('\n');

Ts   = 1/6;
Nsim = 21;

x = zeros(5, Nsim+1);
u = zeros(2, Nsim);

x(:,1) = [15;0;0;0;0];

for k = 1:Nsim
    u(:,k)   = K_stab*x(:,k);                  % or -K_stab*x(:,k) depending on convention
    x(:,k+1) = system_dynamics(x(:,k),u(:,k));
end

g  = x(1,:);
t  = Ts*(0:Nsim);
t_u = Ts*(0:Nsim-1);

figure;
stairs(t_u, u(1,:), 'LineWidth', 1.5); hold on;
stairs(t_u, u(2,:), 'LineWidth', 1.5);
grid on;

xlabel('Time in h','Interpreter','latex');
ylabel('Control inputs','Interpreter','latex');

legend({'$u^{(1)}$','$u^{(2)}$'}, ...
       'Interpreter','latex');

title('Closed-loop control inputs','Interpreter','latex');

cleanfigure;
matlab2tikz('closed_loop_inputs.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');


figure;
plot(t, g, 'LineWidth', 1.5);
grid on;

xlabel('Time in h','Interpreter','latex');
ylabel('Blood glucose level','Interpreter','latex');
title('Closed-loop blood glucose trajectory','Interpreter','latex');

cleanfigure;
matlab2tikz('closed_loop_glucose.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

RES.d.x = x;
RES.d.u = u;


%% --------------------------------------------------------------
% Problem 2(e): Stabilizing gain from noisy state data (dataset 3)
% ---------------------------------------------------------------
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 2(e): Stabilization using noisy state data (dataset 3)\n');
fprintf('--------------------------------------------------------------\n');
fprintf('\n');

% Build Hankel matrices from noisy states and u_off_3
x_off = x_off_3_noisy;     % 5 x (N+1)
u_off = u_off_3;           % 2 x (N+1)

[n, Np1] = size(x_off);
N = Np1 - 1;

x_t = x_off.';             % (N+1) x n
u_t = u_off.';             % (N+1) x m

Hx_n  = H(1, x_t(1:N,   :));   % X_- (noisy)
Hxp_n = H(1, x_t(2:N+1, :));   % X_+ (noisy)
Hu_n  = H(1, u_t(1:N,   :));   % U_- (same inputs)

fprintf('\n--------------------------------------------------------------\n');
fprintf(' Problem 2(e): LMI stabilization using x_off3_noisy}\n');
fprintf('--------------------------------------------------------------\n\n');

% Solver settings (SeDuMi)
sdpopt = sdpsettings('solver','sedumi','verbose',1);

% LMI variables
theta = sdpvar(N,n,'full');
X     = sdpvar(n,n,'symmetric');

% LMI constraints (treat noisy data as if noise-free)
F = [];
F = [F, [X, Hxp_n*theta; theta'*Hxp_n', X] >= 1e-6];
F = [F, X == Hx_n*theta];

st = optimize(F, [], sdpopt);

if st.problem ~= 0
    fprintf('LMI could NOT be certified (solver issue).\n');
    fprintf('YALMIP info: %s\n', st.info);
    K_stabn = [];
else
    ths = value(theta);
    Xs  = value(X);

    if rcond(Xs) < 1e-12
        fprintf('Warning: X is near singular -> K_stabn may be unreliable.\n');
    end

    K_stabn = Hu_n * ths / Xs;

    fprintf('K_stabn = \n');
    disp(K_stabn);

    eig_cl_n = eig(A_true + B_true*K_stabn);
    fprintf('max |eig(A_true + B_true*K_stabn)| = %.4f\n', max(abs(eig_cl_n)));
end


%% --------------------------------------------------------------
% Problem 2(f): Closed-loop simulation with K_stabn
% ---------------------------------------------------------------
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 2(f): Closed-loop simulation with K_{stabn}\n');
fprintf('--------------------------------------------------------------\n');
fprintf('\n');

Ts   = 1/6;
Nsim = 21;

x = zeros(5, Nsim+1);
u = zeros(2, Nsim);

x(:,1) = [15;0;0;0;0];

if isempty(K_stabn)
    error('K_stabn is empty because the LMI solver did not return a valid solution.');
end

for k = 1:Nsim
    u(:,k)   = K_stabn * x(:,k);            % consistent with A + B*K
    x(:,k+1) = system_dynamics(x(:,k),u(:,k));
end

g   = x(1,:);
t   = Ts*(0:Nsim);
t_u = Ts*(0:Nsim-1);

% ---- plot inputs ----
figure;
stairs(t_u, u(1,:), 'LineWidth', 1.5); hold on;
stairs(t_u, u(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time in h','Interpreter','latex');
ylabel('Control inputs','Interpreter','latex');
legend({'$u^{(1)}$','$u^{(2)}$'}, 'Interpreter','latex', 'Location','northeast');
title('Closed-loop control inputs with $K_{\mathrm{stabn}}$','Interpreter','latex');

cleanfigure;
matlab2tikz('closed_loop_inputs_noisy.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% ---- plot glucose ----
figure;
plot(t, g, 'LineWidth', 1.5);
grid on;
xlabel('Time in h','Interpreter','latex');
ylabel('Blood glucose level','Interpreter','latex');
title('Closed-loop blood glucose with $K_{\mathrm{stabn}}$','Interpreter','latex');

cleanfigure;
matlab2tikz('closed_loop_glucose_noisy.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

RES.f.x = x;
RES.f.u = u;


%% --------------------------------------------------------------
% Problem 2(g): Data-based LQR using dataset 3
% --------------------------------------------------------------

fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 2(g): Data-based LQR (dataset 3)\n');
fprintf('--------------------------------------------------------------\n\n');

% weights from the exercise
Q = 2*eye(n);
R = 500*eye(m);

% solver settings
sdpopt = sdpsettings('solver','sedumi','verbose',1);

% ---------------- Step 1: compute P^+ ----------------
P = sdpvar(n,n,'symmetric');

L = Hx_3'*P*Hx_3 - Hxp_3'*P*Hxp_3 - Hx_3'*Q*Hx_3 - Hu_3'*R*Hu_3;

F = [P >= 0, L <= 0];

st = optimize(F, -trace(P), sdpopt);   % maximize trace(P)

if st.problem ~= 0
    error('Step 1 failed (P^+). YALMIP: %s', st.info);
end

Pplus = value(P);

% ---------------- Step 2: compute X_-^\dagger ----------------
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Step 2: compute right inverse X_-^\\dagger\n');
fprintf('--------------------------------------------------------------\n\n');

Ncol = size(Hx_3,2);          % number of columns of X_-
Xdag = sdpvar(Ncol,n,'full'); % right inverse candidate

Lplus = Hx_3'*Pplus*Hx_3 - Hxp_3'*Pplus*Hxp_3 - Hx_3'*Q*Hx_3 - Hu_3'*R*Hu_3;

F = [];
F = [F, Hx_3*Xdag == eye(n)];
F = [F, Lplus*Xdag == zeros(Ncol,n)];

st = optimize(F, [], sdpopt);

% if st.problem ~= 0
%     error('Step 2 failed (X_-^dagger). YALMIP: %s', st.info);
% end

Xdag_s = value(Xdag);

% ---------------- Step 3: K = U_- X_-^\dagger ----------------
K_LQR = Hu_3 * Xdag_s;

fprintf('K_LQR = \n');
disp(K_LQR);

eig_cl = eig(A_true + B_true*K_LQR);
fprintf('max |eig(A_true + B_true*K_LQR)| = %.4f\n', max(abs(eig_cl)));


%% --------------------------------------------------------------
% Problem 2(h): Closed-loop simulation with K_LQR (144 steps)
% ---------------------------------------------------------------
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 2(h): Closed-loop simulation with K_{LQR} (144 steps)\n');
fprintf('--------------------------------------------------------------\n\n');

Ts   = 1/6;     % hours (10 min)
Nsim = 144;

x0 = [15;0;0;0;0];

% preallocate
x = zeros(5, Nsim+1);
u = zeros(2, Nsim);

x(:,1) = x0;

for k = 1:Nsim
    u(:,k)   = K_LQR * x(:,k);                  % u_k = K_LQR x_k
    x(:,k+1) = system_dynamics(x(:,k), u(:,k)); % x_{k+1}
end

% outputs for plotting
g   = x(1,:);              % blood glucose level
t   = Ts*(0:Nsim);         % time for states
t_u = Ts*(0:Nsim-1);       % time for inputs

% ---------------------------------------------------------------
% Plot blood glucose
% ---------------------------------------------------------------
figure;
plot(t, g, 'LineWidth', 1.5);
grid on;
xlabel('Time in h','Interpreter','latex');
ylabel('Blood glucose level','Interpreter','latex');
title('Closed-loop blood glucose trajectory with $K_{\mathrm{LQR}}$','Interpreter','latex');

cleanfigure;
matlab2tikz('closed_loop_glucose_lqr_144.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% ---------------------------------------------------------------
% Plot control inputs
% ---------------------------------------------------------------
figure;
stairs(t_u, u(1,:), 'LineWidth', 1.5); hold on;
stairs(t_u, u(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time in h','Interpreter','latex');
ylabel('Control inputs','Interpreter','latex');
legend({'$u^{(1)}$','$u^{(2)}$'}, 'Interpreter','latex', 'Location','northeast');
title('Closed-loop control inputs with $K_{\mathrm{LQR}}$','Interpreter','latex');

cleanfigure;
matlab2tikz('closed_loop_inputs_lqr_144.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

RES.h.x = x;
RES.h.u = u;

%%
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 2(i): Model-based LQR vs Data-based LQR comparison\n');
fprintf('--------------------------------------------------------------\n\n');

% ---------------------------------------------------------------
% 1) Model-based LQR gain from identified model (IMPORTANT SIGN!)
% MATLAB dlqr returns K such that u = -K x
% Your simulations earlier used u = K x (since you checked A + B*K)
% So we convert: K_LQRmb = -K_dlqr
% ---------------------------------------------------------------
K_dlqr = dlqr(A_id, B_id, Q, R);   % u = -K_dlqr x
K_LQRmb = -K_dlqr;                % u =  K_LQRmb x

fprintf('K_LQRmb = \n'); disp(K_LQRmb);

% quick stability check (with your convention A + B*K)
eig_mb = eig(A_true + B_true*K_LQRmb);
fprintf('max |eig(A_true + B_true*K_LQRmb)| = %.4f\n', max(abs(eig_mb)));

% ---------------------------------------------------------------
% 2) Simulate TRUE system with both controllers
% ---------------------------------------------------------------

% preallocate
x_db = zeros(n, Nsim+1);   u_db = zeros(m, Nsim);   % data-based LQR
x_mb = zeros(n, Nsim+1);   u_mb = zeros(m, Nsim);   % model-based LQR

x_db(:,1) = x0;
x_mb(:,1) = x0;

for k = 1:Nsim
    % data-based LQR
    u_db(:,k)   = K_LQR * x_db(:,k);
    x_db(:,k+1) = system_dynamics(x_db(:,k), u_db(:,k));

    % model-based LQR
    u_mb(:,k)   = K_LQRmb * x_mb(:,k);
    x_mb(:,k+1) = system_dynamics(x_mb(:,k), u_mb(:,k));
end

t  = Ts*(0:Nsim);
t_u = Ts*(0:Nsim-1);

g_db = x_db(1,:);
g_mb = x_mb(1,:);

% ---------------------------------------------------------------
% 3) One plot comparison (blood glucose)
% ---------------------------------------------------------------
figure;
plot(t, g_db, 'LineWidth', 1.5); hold on;
plot(t, g_mb, '--', 'LineWidth', 1.5);
grid on;
xlabel('Time in h','Interpreter','latex');
ylabel('Blood glucose level','Interpreter','latex');
legend({'$g$ with $K_{\mathrm{LQR}}$ (data-based)', ...
        '$g$ with $K_{\mathrm{LQR,mb}}$ (model-based)'}, ...
       'Interpreter','latex','Location','northeast');
title('Data-based vs model-based LQR (true system simulation)','Interpreter','latex');

cleanfigure;
matlab2tikz('compare_glucose_db_vs_mb.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% ---------------------------------------------------------------
% Optional: compare inputs too
% ---------------------------------------------------------------
figure;
stairs(t_u, u_db(1,:), 'LineWidth', 1.2); hold on;
stairs(t_u, u_mb(1,:), '--', 'LineWidth', 1.2);
stairs(t_u, u_db(2,:), 'LineWidth', 1.2);
stairs(t_u, u_mb(2,:), '--', 'LineWidth', 1.2);
grid on;
xlabel('Time in h','Interpreter','latex');
ylabel('Control inputs','Interpreter','latex');
legend({'$u^{(1)}$ data-based', '$u^{(1)}$ model-based', ...
        '$u^{(2)}$ data-based', '$u^{(2)}$ model-based'}, ...
       'Interpreter','latex','Location','northeast');
title('Control inputs: data-based vs model-based LQR','Interpreter','latex');

cleanfigure;
matlab2tikz('compare_inputs_db_vs_mb.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

%% --------------------------------------------------------------
% Problem 2(j): Data-based LQR with noisy states (as if noise-free)
% --------------------------------------------------------------

fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 2(j): Data-based LQR with noisy states (as if noise-free)\n');
fprintf('--------------------------------------------------------------\n\n');

% ---------------------------------------------------------------
% Build X_-, X_+, U_- from (u_off_3, x_off_3_noisy)
% NOTE: we use L = 1 Hankel, so it is just stacking the samples
% ---------------------------------------------------------------
x_tn = x_off_3_noisy.';   % (N+1) x n
u_t  = u_off_3.';         % (N+1) x m

Nn = size(x_tn,1) - 1;    % N

Hx_n  = H(1, x_tn(1:Nn,   :));   % X_- (noisy)
Hxp_n = H(1, x_tn(2:Nn+1, :));   % X_+ (noisy)
Hu_n  = H(1, u_t(1:Nn,    :));   % U_- (same u)

% solver
sdpopt = sdpsettings('solver','sedumi','verbose',1);

% ---------------------------------------------------------------
% Step 1 (slides): compute P^+  via L(P) <= 0
% ---------------------------------------------------------------
P = sdpvar(n,n,'symmetric');

L = Hx_n'*P*Hx_n - Hxp_n'*P*Hxp_n - Hx_n'*Q*Hx_n - Hu_n'*R*Hu_n;

F = [P >= 0, L <= 0];

st = optimize(F, -trace(P), sdpopt);

if st.problem ~= 0
    error('2(j) Step 1 failed (P^+). YALMIP: %s', st.info);
end

Pplus = value(P);

% ---------------------------------------------------------------
% Step 2 (slides): compute right inverse X_-^\dagger
% ---------------------------------------------------------------
Ncol = size(Hx_n,2);
Xdag = sdpvar(Ncol,n,'full');

Lplus = Hx_n'*Pplus*Hx_n - Hxp_n'*Pplus*Hxp_n - Hx_n'*Q*Hx_n - Hu_n'*R*Hu_n;

F = [];
F = [F, Hx_n*Xdag == eye(n)];
F = [F, Lplus*Xdag == zeros(Ncol,n)];

st = optimize(F, [], sdpopt);

Xdag_s = value(Xdag);

% controller gain
K_LQRn = Hu_n * Xdag_s;

fprintf('K_LQRn = \n');
disp(K_LQRn);

% optional stability check with TRUE matrices (your sign convention A + B*K)
eig_cl = eig(A_true + B_true*K_LQRn);
fprintf('max |eig(A_true + B_true*K_LQRn)| = %.4f\n', max(abs(eig_cl)));

%% --------------------------------------------------------------
% Problem 2(k): Closed-loop simulation with K_LQRn (144 steps)
% ---------------------------------------------------------------
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 2(k): Closed-loop simulation with K_{LQRn} (144 steps)\n');
fprintf('--------------------------------------------------------------\n\n');

Ts   = 1/6;     % hours (10 min)
Nsim = 144;

x0 = [15;0;0;0;0];

% preallocate
x = zeros(n, Nsim+1);
u = zeros(m, Nsim);

x(:,1) = x0;

% ---- simulate TRUE system ----
for k = 1:Nsim
    u(:,k)   = K_LQRn * x(:,k);                  % u_k = K_LQRn x_k
    x(:,k+1) = system_dynamics(x(:,k), u(:,k));  % TRUE dynamics
end

% time vectors
t   = Ts*(0:Nsim);
t_u = Ts*(0:Nsim-1);

% output of interest
g = x(1,:);   % glucose

% ---- stability check (with your convention A + B*K) ----
eig_cl_n = eig(A_true + B_true*K_LQRn);
fprintf('eig(A_true + B_true*K_LQRn) = \n'); disp(eig_cl_n);
fprintf('max |eig| = %.6f\n', max(abs(eig_cl_n)));

% ---- plot glucose ----
figure;
plot(t, g, 'LineWidth', 1.5);
grid on;
xlabel('Time in h','Interpreter','latex');
ylabel('Blood glucose level','Interpreter','latex');
title('Closed-loop blood glucose with $K_{\mathrm{LQRn}}$','Interpreter','latex');

cleanfigure;
matlab2tikz('closed_loop_glucose_lqrn_144.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% ---- plot inputs ----
figure;
stairs(t_u, u(1,:), 'LineWidth', 1.5); hold on;
stairs(t_u, u(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time in h','Interpreter','latex');
ylabel('Control inputs','Interpreter','latex');
legend({'$u^{(1)}$','$u^{(2)}$'}, 'Interpreter','latex', 'Location','northeast');
title('Closed-loop inputs with $K_{\mathrm{LQRn}}$','Interpreter','latex');

cleanfigure;
matlab2tikz('closed_loop_inputs_lqrn_144.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');
RES.k.x = x;
RES.k.u = u;


%% --------------------------------------------------------------
% Comparison of control inputs and glucose (Problem 2(l))
% --------------------------------------------------------------
fprintf('\n==============================================================\n');
fprintf(' Problem 2(l): Comparison of control inputs and glucose drift\n');
fprintf('==============================================================\n\n');

labels = {'(d) K_{stab}', '(f) K_{stabn}', '(h) K_{LQR}', '(k) K_{LQRn}'};
keys   = {'d','f','h','k'};

for ii = 1:length(keys)
    key = keys{ii};

    if ~isfield(RES, key)
        continue
    end

    x = RES.(key).x;
    u = RES.(key).u;

    g = x(1,:);
    du = diff(u,1,2);

    fprintf('%s\n', labels{ii});
    fprintf('  g(1) = %.6g, g(end) = %.6g, g(end)-g(1) = %.6g\n', ...
        g(1), g(end), g(end)-g(1));

    for j = 1:size(u,1)
        fprintf('  max |u_%d| = %.6g,   max |Delta u_%d| = %.6g\n', ...
            j, max(abs(u(j,:))), j, max(abs(du(j,:))));
    end
    fprintf('\n');
end

fprintf('==============================================================\n\n');


%% --------------------------------------------------------------
% Problem 2(m): Only insulin controlled, carbs fixed to u_sim from Problem 1
% --------------------------------------------------------------
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 2(m): Insulin-only feedback + regular eating (u_sim)\n');
fprintf('--------------------------------------------------------------\n\n');

% --- load u_sim from problem 1 ---
S1 = load('data_problem1.mat');   % contains u_sim
u_sim = S1.u_sim;                 % should be 2 x L (or L x 2)
clear S1;

% make sure u_sim is 2 x Nsim
if size(u_sim,1) ~= 2 && size(u_sim,2) == 2
    u_sim = u_sim.';              % convert to 2 x L if needed
end

Ts   = 1/6;      % hours
Nsim = 144;
n = 5; m = 2;

% we need the second input component over 144 steps
% --- extend u_sim periodically if it is too short ---
L = size(u_sim,2);

if L < Nsim
    reps = ceil(Nsim / L);
    u_sim_ext = repmat(u_sim, 1, reps);
else
    u_sim_ext = u_sim;
end

u_carbs = u_sim_ext(2,1:Nsim);   % regular eating pattern


% --- build Kstab2: only insulin is feedback-controlled ---
Kstab2 = K_stab;
Kstab2(2,:) = 0;                  % second input not controlled by feedback

fprintf('Kstab2 = \n'); disp(Kstab2);

% --- simulate TRUE system ---
x = zeros(n, Nsim+1);
u = zeros(m, Nsim);

x(:,1) = [15;0;0;0;0];

for k = 1:Nsim
    u_fb = Kstab2 * x(:,k);        % u = Kx, but 2nd row is zero
    u(1,k) = u_fb(1);              % insulin input (feedback)
    u(2,k) = u_carbs(k);           % carbs input (preset regular eating)
    x(:,k+1) = system_dynamics(x(:,k), u(:,k));
end

g   = x(1,:);
t   = Ts*(0:Nsim);
t_u = Ts*(0:Nsim-1);

% --- plot glucose ---
figure;
plot(t, g, 'LineWidth', 1.5);
grid on;
xlabel('Time in h','Interpreter','latex');
ylabel('Blood glucose level','Interpreter','latex');
title('Glucose with insulin-only feedback ($K_{\mathrm{stab2}}$)','Interpreter','latex');

cleanfigure;
matlab2tikz('problem2m_glucose.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% --- plot inputs ---
figure;
stairs(t_u, u(1,:), 'LineWidth', 1.5); hold on;
stairs(t_u, u(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time in h','Interpreter','latex');
ylabel('Control inputs','Interpreter','latex');
legend({'$u^{(1)}$ insulin (feedback)', '$u^{(2)}$ carbs (preset $u_{\mathrm{sim}}$)'}, ...
       'Interpreter','latex','Location','northeast');
title('Inputs with insulin-only feedback','Interpreter','latex');

cleanfigure;
matlab2tikz('problem2m_inputs.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% --- stability check (IMPORTANT NOTE BELOW) ---
eig_cl2 = eig(A_true + B_true*Kstab2);
fprintf('eig(A_true + B_true*Kstab2) = \n'); disp(eig_cl2);
fprintf('max |eig| = %.12f\n', max(abs(eig_cl2)));

% optional drift indicators
fprintf('g(1)=%.6g, g(end)=%.6g, g(end)-g(1)=%.6g\n', g(1), g(end), g(end)-g(1));
fprintf('max |u1|=%.6g, max |u2|=%.6g\n', max(abs(u(1,:))), max(abs(u(2,:))));
