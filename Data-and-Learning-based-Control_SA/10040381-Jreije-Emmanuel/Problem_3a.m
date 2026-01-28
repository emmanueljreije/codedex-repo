clear; clc; close all;
% addpath(genpath('C:\Users\CompuTop\seadrive_root\Emmanuel\Meine Bibliotheken\My Library\MATLAB-Libs\sedumi-master'))
addpath(genpath('C:\Users\CompuTop\seadrive_root\Emmanuel\Meine Bibliotheken\My Library\MATLAB-Libs\YALMIP-master'))
addpath(genpath("C:\Users\CompuTop\seadrive_root\Emmanuel\Meine Bibliotheken\My Library\MATLAB-Libs\sedumi"))
addpath(genpath("C:\Users\CompuTop\seadrive_root\Emmanuel\Meine Bibliotheken\My Library\MATLAB-Libs\matlab2tikz-master"))

%% --------------------------------------------------------------
% Problem 3(a): Pendulum setpoint stabilization (full model knowledge)
% ---------------------------------------------------------------
fprintf('\n');
fprintf('--------------------------------------------------------------\n');
fprintf(' Problem 3(a): Pendulum setpoint stabilization (Euler)\n');
fprintf('--------------------------------------------------------------\n\n');

% constants
g = 9.81;    % m/s^2
l = 10;      % m

% controller gains from pole placement: (lambda+6)^2 = lambda^2 + 12 lambda + 36
k1 = 36;
k2 = 12;

% simulation settings (given)
Ts   = 0.01;    % s
Nsim = 2000;

x0 = [1; 0];

% preallocate
x = zeros(2, Nsim+1);
u = zeros(1, Nsim);

x(:,1) = x0;

for k = 1:Nsim
    x1 = x(1,k);
    x2 = x(2,k);

    % control law (full model knowledge)
    u(k) = (g/l)*sin(x1) - k1*x1 - k2*x2;

    % true dynamics
    xdot1 = x2;
    xdot2 = -(g/l)*sin(x1) + u(k);

    % Euler step
    x(:,k+1) = x(:,k) + Ts*[xdot1; xdot2];
end

% time vectors
t   = Ts*(0:Nsim);      % time for states
t_u = Ts*(0:Nsim-1);    % time for inputs

% ---------------------------------------------------------------
% Plot states
% ---------------------------------------------------------------
figure;
plot(t, x(1,:), 'LineWidth', 1.5); hold on;
plot(t, x(2,:), 'LineWidth', 1.5);
grid on;
xlabel('Time in s','Interpreter','latex');
ylabel('States','Interpreter','latex');
legend({'$x_1$ (angle)',' $x_2$ (angular velocity)'}, ...
    'Interpreter','latex', 'Location','northeast');
title('Pendulum closed-loop states (Problem 3a)','Interpreter','latex');

cleanfigure;
matlab2tikz('p3a_states.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% ---------------------------------------------------------------
% Plot control input
% ---------------------------------------------------------------
figure;
plot(t_u, u, 'LineWidth', 1.5);
grid on;
xlabel('Time in s','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
title('Pendulum control input (Problem 3a)','Interpreter','latex');

cleanfigure;
matlab2tikz('p3a_input.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

% store results
RES.p3a.x = x;
RES.p3a.u = u;
RES.p3a.t = t;