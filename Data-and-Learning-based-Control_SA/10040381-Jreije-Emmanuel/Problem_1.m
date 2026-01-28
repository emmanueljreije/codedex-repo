clear; clc; close all;
addpath(genpath('matlab2tikz-master'))
load('data_problem1.mat');   

% ------------------------------------------------------------------
% Part a) Plots of the control inputs
% ------------------------------------------------------------------
Ts = 10;    % sampling time in minutes

% time vector
N = size(u_off_1, 2);  
t = (0:N-1) * Ts;       

font_title = 14;
font_label = 14;
font_legend = 10;
font_ticks = 10;
font_sgtitle = 12;

c1 = [0.1 0.2 0.8];   % blue
c2 = [0.9 0.4 0.1];   % orange
c3 = [0.1 0.6 0.2];   % green

figure;
tl = tiledlayout(6,1);          % 6 rows, 1 column
tl.TileSpacing = 'tight';     % or 'loose' for more gaps
tl.Padding     = 'none';     % less empty border

%% 1) u_off_1, input 1
nexttile;
plot(t, u_off_1(1,:), 'LineWidth', 1, 'Color', c1);
ylabel('$u_{\mathrm{off},1}^{(1)}$','Interpreter','latex','FontSize',font_label);
set(gca,'FontSize',font_ticks);
grid on;
axis tight;

%% 2) u_off_1, input 2
nexttile;
plot(t, u_off_1(2,:), 'LineWidth', 1, "Color", c1);
ylabel('$u_{\mathrm{off},1}^{(2)}$','Interpreter','latex','FontSize',font_label);
set(gca,'FontSize',font_ticks);
grid on;
axis tight;

%% 3) u_off_2, input 1
nexttile;
plot(t, u_off_2(1,:), 'LineWidth', 1, 'color', c2);
ylabel('$u_{\mathrm{off},2}^{(1)}$','Interpreter','latex','FontSize',font_label);
set(gca,'FontSize',font_ticks);
grid on;
axis tight;

%% 4) u_off_2, input 2
nexttile;
plot(t, u_off_2(2,:), 'LineWidth', 1, 'color', c2);
ylabel('$u_{\mathrm{off},2}^{(2)}$','Interpreter','latex','FontSize',font_label);
set(gca,'FontSize',font_ticks);
grid on;
axis tight;

%% 5) u_off_3, input 1
nexttile;
plot(t, u_off_3(1,:), 'LineWidth', 1, 'Color', c3);
ylabel('$u_{\mathrm{off},3}^{(1)}$','Interpreter','latex','FontSize',font_label);
set(gca,'FontSize',font_ticks);
grid on;
axis tight;

%% 6) u_off_3, input 2  (only here we add x-label)
nexttile;
plot(t, u_off_3(2,:), 'LineWidth', 1, 'Color', c3);
xlabel('Time in min','Interpreter','latex','FontSize',font_label);
ylabel('$u_{\mathrm{off},3}^{(2)}$','Interpreter','latex','FontSize',font_label);
set(gca,'FontSize',font_ticks);
grid on;
axis tight;

sgtitle('Offline control input sequences', ...
        'Interpreter','latex','FontSize',font_sgtitle);

cleanfigure;
matlab2tikz('offline_inputs.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');


% ------------------------------------------------------------------
% Part b) Checking for persistently exiting control input trajectory
% ------------------------------------------------------------------

L = 72;          % number of simulated time steps (given)
n = 5;           % number of states: [g I1 I2 D1 D2]
L_PE = L + n;    % order required by Willems' lemma

datasets = {u_off_1, u_off_2, u_off_3};
names    = {'u_{off,1}', 'u_{off,2}', 'u_{off,3}'};

for k = 1:numel(datasets)
    u = datasets{k};   % u is 2 × N (2 inputs: insulin + carbs)

    [isPE, rH, nRows] = is_persistently_exciting(u, L_PE);

    if isPE
        fprintf('The control input %s IS persistently exciting of order L+n = %d (rank(H) = %d, rows = %d).\n', ...
                names{k}, L_PE, rH, nRows);
    else
        fprintf('The control input %s is NOT persistently exciting of order L+n = %d (rank(H) = %d, rows = %d).\n', ...
                names{k}, L_PE, rH, nRows);
    end
end


% ------------------------------------------------------------------
% Part c) Simulation via Willems' fundamental lemma
% ------------------------------------------------------------------

% dimensions
m = 2;       % number of inputs (insulin, carbs)
p = 1;       % number of outputs (blood glucose)
n = 5;       % number of states [g I1 I2 D1 D2]
L = 72;      % number of future time steps to simulate
Tini = n;    % length of initial trajectory (must be ≥ system order)

% u_off_3 IS persistently exciting of order L+n = 77
u_off = u_off_3;    % 2 x N_off
y_off = y_off_3;    % 1 x N_off

fprintf('Using dataset u_{off,3} / y_{off,3} as PE data for Willems simulation.\n');

% ------------------------------------------------------------------
% 1) Build Hankel matrices for (u_off, y_off)
% ------------------------------------------------------------------
N_off   = size(u_off, 2);       % length of offline data
L_total = Tini + L;             % total window length

if N_off < L_total
    error('Offline dataset u_{off,3} is too short: need at least Tini + L samples.');
end

% H expects data as N x m (time along rows)
Hu = H(L_total, u_off.');   % size: (m*L_total) x (N_off - L_total + 1)
Hy = H(L_total, y_off.');   % size: (p*L_total) x (N_off - L_total + 1)

% split past/future blocks
U_p = Hu(1:m*Tini, :);             % past inputs (Tini steps)
U_f = Hu(m*Tini+1:end, :);         % future inputs (L steps)

Y_p = Hy(1:p*Tini, :);             % past outputs (Tini steps)
Y_f = Hy(p*Tini+1:end, :);         % future outputs (L steps)

% ------------------------------------------------------------------
% 2) Build desired trajectory pieces (initial 0, given future input u_sim)
% ------------------------------------------------------------------

% initial input/output trajectory at zero
u_ini = zeros(m*Tini, 1);          % stacked [u_0; ...; u_{Tini-1}]
y_ini = zeros(p*Tini, 1);          % stacked [y_0; ...; y_{Tini-1}]

% online input trajectory for future (use first L samples of u_sim)
u_f_mat = u_sim(:, 1:L);           % 2 x L
u_f     = u_f_mat(:);              % vectorized (mL x 1)

% ------------------------------------------------------------------
% 3) Solve for g using Willems' lemma:
%    [U_p; Y_p; U_f] * g = [u_ini; y_ini; u_f]
% ------------------------------------------------------------------
M   = [U_p; Y_p; U_f];             % big data matrix
rhs = [u_ini; y_ini; u_f];         % desired trajectory

% least-squares solution
g = M \ rhs;

% future outputs from data
y_f_stack = Y_f * g;               % size pL x 1
y_hat     = reshape(y_f_stack, p, L);   % 1 x L

% ground truth simulation trajectory (given in data_problem1.mat)
y_true = y_sim(1:L);               % 1 x L

% ------------------------------------------------------------------
% 4) Plot comparison in hours
% ------------------------------------------------------------------
t_hours = (0:L-1) * Ts / 60;       % convert minutes to hours

figure;
plot(t_hours, y_true, 'LineWidth', 1.5); hold on;
plot(t_hours, y_hat, '--', 'LineWidth', 1.5);
grid on;

xlabel('Time in h','Interpreter','latex','FontSize',font_label);
ylabel('Blood glucose level','Interpreter','latex','FontSize',font_label);

legend({'$y_{\mathrm{sim}}$', '$\hat{y}_{\mathrm{Willems}}$'}, ...
       'Interpreter','latex','FontSize',font_legend, 'Location','northwest');

set(gca,'FontSize',font_ticks);

sgtitle('Blood glucose simulation via Willems'' fundamental lemma', ...
        'Interpreter','latex','FontSize',font_sgtitle);

% optional tikz export
cleanfigure;
matlab2tikz('willems_simulation.tex', ...
    'standalone', false, ...
    'showInfo', false, ...
    'parseStrings', false, ...
    'floatFormat', '%.6f', ...
    'height', '\figureheight', ...
    'width',  '\figurewidth');

