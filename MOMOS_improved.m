% The MOMOS model simulation with interactive case selection
clear all; close all; clc;
tic;  % Start execution timer

%% PDE System Parameters
Du = 0.6;
Dv = 0.6;
k1 = 0.4;
k2 = 0.6;
c  = 0.8;

%% Ask user to select one of the 4 cases
disp('Select one of the following cases:');
disp('1 - Case 1: q = 0.0433,  beta = 0.806');
disp('2 - Case 2: q = 0.061122, beta = 1.01668');
disp('3 - Case 3: q = 0.0196639, beta = 0.474095');
disp('4 - Case 4: q = 0.0804361, beta = 1.23535');
case_choice = input('Enter the number of the case you want to run (1-4): ');

% Assign parameters based on user selection
switch case_choice
    case 1
        q = 0.0433; beta = 0.806; sigma = 10;
    case 2
        q = 0.061122; beta = 1.01668; sigma = 12;
    case 3
        q = 0.0196639; beta = 0.474095; sigma = 12;
    case 4
        q = 0.0804361; beta = 1.23535; sigma = 7;
    otherwise
        error('Invalid selection. Please choose a number between 1 and 4.');
end

%% Define reaction terms
f = @(u, v) -k1*u - q*u.*abs(u) + k2*v;
g = @(u, v)  k1*u - k2*v + c;

%% Domain setup
L = 15;
h = 0.2;
n = round(L/h) + 1;
x = linspace(0, L, n);
y = linspace(0, L, n);
[xx, yy] = meshgrid(x, y);

%% Time discretization
t0 = 0; tf = 500; ht = 0.01;
t = t0:ht:tf;
nt = length(t);

%% Laplacian (2D) with periodic BCs
N = n^2;
I = speye(n); evec = ones(n,1);
Lap = spdiags([evec, -2*evec, evec], [1, 0, -1], n, n);
Lap(1,end) = 1; Lap(end,1) = 1;
A = ((1/h)^2) * kron(I, Lap) + ((1/h)^2) * kron(Lap, I);

%% Chemotaxis discretization
A1 = spdiags([evec evec], [0 1], n, n);  A1(end,1) = 1;
A2 = spdiags([evec evec], [0 -1], n, n); A2(1,end) = 1;
B1 = spdiags([-evec evec], [0 1], n, n); B1(end,1) = 1;
B2 = spdiags([evec -evec], [0 -1], n, n); B2(1,end) = -1;

A1 = 0.5 * A1; A2 = 0.5 * A2;

A1x = (1/h) * kron(I, A1);   A1y = (1/h) * kron(A1, I);
A2x = (1/h) * kron(I, A2);   A2y = (1/h) * kron(A2, I);
B1x = (1/h) * kron(I, B1);   B1y = (1/h) * kron(B1, I);
B2x = (1/h) * kron(I, B2);   B2y = (1/h) * kron(B2, I);

%% Initial Conditions
rng(42);  % for reproducibility

ue = sqrt(c / q);
ve = (k1 / k2) * sqrt(c / q) + c / k2;

switch case_choice
    case {1, 3}
        % Random perturbation around equilibrium
        U0 = ue + sigma * rand(n, n);
        V0 = ve + sigma * rand(n, n);
    case {2, 4}
        % Band with high constant value in top 30%
        U0 = ue * (1 + sigma * rand(n, n));
        thickness = ceil(0.3 * n);
        constant_value = 40;
        U0(1:thickness, :) = constant_value;
        V0 = ve * ones(n, n);
end

% Flatten for vector operations
u = U0(:);
v = V0(:);

%% Preprocessing for implicit scheme
In = speye(N);
Mu = In - Du * ht * A;
Mv = (1 + k2 * ht) * In - Dv * ht * A;
[Lmu, Umu] = lu(Mu);
[Lmv, Umv] = lu(Mv);
ht_beta = ht * beta;

%% Time-stepping loop
for k = 1 : nt - 1
    % Update v
    v = Umv \ (Lmv \ (v + ht * (k1 * u + c)));

    % Chemotaxis term
    chemotaxis = ...
        (A1x * u) .* (B1x * v) - (A2x * u) .* (B2x * v) + ...
        (A1y * u) .* (B1y * v) - (A2y * u) .* (B2y * v);

    % Update u
    u = Umu \ (Lmu \ (u + ht * f(u, v) - ht_beta * chemotaxis));
end

%% Report execution time
fprintf('\nExecution Time: %.2f seconds\n', toc);

%% Plot final results
figure;
pcolor(xx, yy, reshape(u, n, n));
shading interp; colormap('jet'); colorbar;
xlabel('x', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 14, 'FontWeight', 'bold');
title('u at T = 500', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');

figure;
pcolor(xx, yy, reshape(v, n, n));
shading interp; colormap('jet'); colorbar;
xlabel('x', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 14, 'FontWeight', 'bold');
title('v at T = 500', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
