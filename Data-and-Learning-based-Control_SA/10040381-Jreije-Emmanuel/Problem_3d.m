%% --------------------------------------------------------------
% Problem 3(d): GP-MRAC tracking Van der Pol reference (Euler)
% ---------------------------------------------------------------
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 3(d): GP-MRAC tracks Van der Pol reference (Euler)\n');
fprintf('--------------------------------------------------------------\n\n');

% constants
g = 9.81;    % m/s^2
l = 10;      % m

% feedback gains from (a)
k1 = 36;
k2 = 12;

% Van der Pol parameter
mu = 1.5;

% GP kernel hyperparameter 
ell = 1;                 % length-scale l=1
sigma_n2 = 1e-6;         % nugget for numerical stability
Mmax = 10;               % use only last 10 data pairs 

% simulation settings 
Ts   = 0.01;             % s
Nsim = 2000;

% initial conditions 
x0  = [1; 0];            % true system
xr0 = [1.5; 0];          % reference system

% preallocate
x  = zeros(2, Nsim+1);
xr = zeros(2, Nsim+1);
u  = zeros(1, Nsim);

% store mismatch info (useful for report/debug)
Delta_true = zeros(1, Nsim);   % computed training targets
Delta_hat  = zeros(1, Nsim);   % GP mean used in control

x(:,1)  = x0;
xr(:,1) = xr0;

% GP training buffers (store last Mmax points)
Xtr = zeros(2, 0);   % each column is a 2D state x
Ytr = zeros(1, 0);   % mismatch targets

for k = 1:Nsim
    % current states
    xk  = x(:,k);
    xrk = xr(:,k);

    x1  = xk(1);   x2  = xk(2);
    xr1 = xrk(1);  xr2 = xrk(2);

    % ------------------------------
    % reference dynamics (Van der Pol)
    % ------------------------------
    xr_dot1 = xr2;
    xr_dot2 = mu*(1 - xr1^2)*xr2 - xr1;

    % tracking error
    e  = xk - xrk;
    e1 = e(1);
    e2 = e(2);

    % ------------------------------
    % GP prediction: mu_Delta(xk)
    % ------------------------------
    if isempty(Ytr)
        muDelta = 0;   % initialize v_a = 0 
    else
        M = size(Xtr,2);

        % K(Xtr,Xtr)
        K = zeros(M,M);
        for i = 1:M
            for j = 1:M
                d = Xtr(:,i) - Xtr(:,j);
                K(i,j) = exp(-(d.'*d)/(2*ell^2));
            end
        end
        K = K + sigma_n2*eye(M);

        % k_*(Xtr,xk)
        kstar = zeros(M,1);
        for i = 1:M
            d = Xtr(:,i) - xk;
            kstar(i) = exp(-(d.'*d)/(2*ell^2));
        end

        alpha = K \ (Ytr.');       % solve K*alpha = Y
        muDelta = kstar.' * alpha; % scalar
    end
    Delta_hat(k) = muDelta;

    % ------------------------------
    % estimated model 
    % x2dot_hat = fhat(x) + u, with ghat=1
    % fhat(x) = -(g/l) sin(x1) + 30 cos(x1)
    % ------------------------------
    fhat = -(g/l)*sin(x1) + 30*cos(x1);

    % desired x2dot (tracking)
    x2dot_des = xr_dot2 - k1*e1 - k2*e2;

    % control law with GP compensation
    % u = x2dot_des - fhat - muDelta
    u(k) = x2dot_des - fhat - muDelta;

    % ------------------------------
    % true system dynamics 
    % ------------------------------
    xdot1 = x2;
    xdot2 = -(g/l)*sin(x1) + u(k);

    % Euler step: plant
    x(:,k+1) = x(:,k) + Ts*[xdot1; xdot2];

    % Euler step: reference
    xr(:,k+1) = xr(:,k) + Ts*[xr_dot1; xr_dot2];

    % ------------------------------
    % build GP training target Delta
    % ------------------------------
    x2dot_true_meas = (x(2,k+1) - x(2,k)) / Ts;  % finite difference
    x2dot_hat = fhat + u(k);                     % model prediction
    delta_k = x2dot_true_meas - x2dot_hat;       % mismatch
    Delta_true(k) = delta_k;

    % append to training data (keep last Mmax)
    Xtr = [Xtr, xk];
    Ytr = [Ytr, delta_k];

    if size(Xtr,2) > Mmax
        Xtr = Xtr(:, end-Mmax+1:end);
        Ytr = Ytr(:, end-Mmax+1:end);
    end
end

% time vectors
t   = Ts*(0:Nsim);
t_u = Ts*(0:Nsim-1);

% ---------------------------------------------------------------
% Plot: x1 vs xr1
% ---------------------------------------------------------------
figure;
plot(t, x(1,:), 'LineWidth', 1.5); hold on;
plot(t, xr(1,:), 'LineWidth', 1.5);
grid on;
xlabel('Time in s','Interpreter','latex');
ylabel('$x_1$','Interpreter','latex');
legend({'$x_1$ (pendulum)',' $x_{r1}$ (reference)'}, ...
    'Interpreter','latex', 'Location','northeast');
title('Problem 3(d): Tracking $x_{r1}$ with GP-MRAC','Interpreter','latex');

cleanfigure;
matlab2tikz('p3d_x1_tracking.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% ---------------------------------------------------------------
% Plot: x2 vs xr2
% ---------------------------------------------------------------
figure;
plot(t, x(2,:), 'LineWidth', 1.5); hold on;
plot(t, xr(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time in s','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
legend({'$x_2$ (pendulum)',' $x_{r2}$ (reference)'}, ...
    'Interpreter','latex', 'Location','northeast');
title('Problem 3(d): Tracking $x_{r2}$ with GP-MRAC','Interpreter','latex');

cleanfigure;
matlab2tikz('p3d_x2_tracking.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% ---------------------------------------------------------------
% Plot: control input
% ---------------------------------------------------------------
figure;
plot(t_u, u, 'LineWidth', 1.5);
grid on;
xlabel('Time in s','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
title('Problem 3(d): Control input','Interpreter','latex');

cleanfigure;
matlab2tikz('p3d_input.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% ---------------------------------------------------------------
% Plot: mismatch learning (optional but useful)
% ---------------------------------------------------------------
figure;
plot(t_u, Delta_true, 'LineWidth', 1.5); hold on;
plot(t_u, Delta_hat,  'LineWidth', 1.5);
grid on;
xlabel('Time in s','Interpreter','latex');
ylabel('$\Delta$','Interpreter','latex');
legend({'$\Delta$ (measured)',' $\mu_\Delta(x)$ (GP mean)'}, ...
    'Interpreter','latex', 'Location','northeast');
title('Problem 3(d): GP mismatch learning','Interpreter','latex');

cleanfigure;
matlab2tikz('p3d_delta_learning.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% store results
RES.p3d.x = x;
RES.p3d.xr = xr;
RES.p3d.u = u;
RES.p3d.Delta_true = Delta_true;
RES.p3d.Delta_hat  = Delta_hat;
RES.p3d.t = t;
