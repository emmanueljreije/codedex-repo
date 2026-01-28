%% --------------------------------------------------------------
% Problem 3(b): Tracking a Van der Pol reference (full model knowledge)
% ---------------------------------------------------------------
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 3(b): Pendulum tracks Van der Pol reference (Euler)\n');
fprintf('--------------------------------------------------------------\n\n');

% constants
g = 9.81;    % m/s^2
l = 10;      % m

% same controller gains as in (a)
k1 = 36;
k2 = 12;

% Van der Pol parameter
mu = 1.5;

% simulation settings
Ts   = 0.01;   % s
Nsim = 2000;

% initial conditions
x0  = [1; 0];
xr0 = [1.5; 0];

% preallocate
x  = zeros(2, Nsim+1);
xr = zeros(2, Nsim+1);
u  = zeros(1, Nsim);

x(:,1)  = x0;
xr(:,1) = xr0;

for k = 1:Nsim
    % current states
    x1  = x(1,k);   x2  = x(2,k);
    xr1 = xr(1,k);  xr2 = xr(2,k);

    % reference dynamics (Van der Pol)
    xr_dot1 = xr2;
    xr_dot2 = mu*(1 - xr1^2)*xr2 - xr1;   % = d/dt (xr2)

    % tracking error
    e1 = x1 - xr1;
    e2 = x2 - xr2;

    % tracking control law
    % u = (g/l) sin(x1) + v,  v = xr_dot2 - k1 e1 - k2 e2
    u(k) = (g/l)*sin(x1) + (xr_dot2 - k1*e1 - k2*e2);

    % true pendulum dynamics
    xdot1 = x2;
    xdot2 = -(g/l)*sin(x1) + u(k);

    % Euler step: plant
    x(:,k+1) = x(:,k) + Ts*[xdot1; xdot2];

    % Euler step: reference
    xr(:,k+1) = xr(:,k) + Ts*[xr_dot1; xr_dot2];
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
title('Tracking of $x_{r1}$ (Problem 3b)','Interpreter','latex');

cleanfigure;
matlab2tikz('p3b_x1_tracking.tex', ...
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
title('Tracking of $x_{r2}$ (Problem 3b)','Interpreter','latex');

cleanfigure;
matlab2tikz('p3b_x2_tracking.tex', ...
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
title('Control input for tracking (Problem 3b)','Interpreter','latex');

cleanfigure;
matlab2tikz('p3b_input.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% store results
RES.p3b.x  = x;
RES.p3b.xr = xr;
RES.p3b.u  = u;
RES.p3b.t  = t;
