%% -----------------------------------------------------------
%  Skript: solve_fuel_minimization_corrected.m
%  Zadanie F: Optimálne riadenie spotreby paliva
%             pre N = 30, s danými A, b, x_fin.
%
%  V tomto skripte sú OPRAVENÉ znamienka v nerovnostiach
%  pre z_t, aby správne platilo z_t >= f(u(t)).
%% -----------------------------------------------------------

clear; clc;

%% (1) PARAMETRE ÚLOHY
N = 30;

A_mat = [ -1   0.4   0.8
           1   0     0
           0   1     0 ];
b_vec = [1; 0; 0.3];
x_fin = [7; 2; -6];

n = size(A_mat,1);  % Rozmer stavu (tu n=3).

%% (2) VYTVORENIE MATICE H PRE x(N) = H * u
%    "u" zatiaľ len formálne ako u(0..N-1).
%    Neskôr zohľadníme, že u(t)=u_t^+ - u_t^-.
H = zeros(n, N);
for k = 0:(N-1)
    H(:,k+1) = (A_mat^(N-1-k)) * b_vec; 
end

%% (3) FORMULÁCIA LP S PRIDAVNÝMI PREMENNÝMI
%  Premenné:
%    pre každý t=0..N-1:
%       u_t^+ >= 0,  u_t^- >= 0,  pričom u(t)=u_t^+ - u_t^- 
%       z_t >= 0 (pomocná premenná na f(u(t)) = max(|u|, 2|u| - 1))
%
%  Celkový vektor X má dĺžku 3*N:
%   [u_0^+, u_0^-, z_0,  u_1^+, u_1^-, z_1, ..., u_{N-1}^+, u_{N-1}^-, z_{N-1}].
%
%  Cieľ: min ∑ z_t.

% (3a) Definícia vektora c (rozmer 3*N), 
%      kde c má na pozíciách z_t hodnotu 1, inak 0.
c = zeros(3*N, 1);
for t = 0:(N-1)
    idx_z = 3*t + 3;
    c(idx_z) = 1;
end

% (3b) Rovnica x(N)=x_fin => H*u = x_fin,
%      ale u(t)=u_t^+ - u_t^- => potrebujeme maticu M (N x 3N), 
%      aby u = M*X.
M = zeros(N, 3*N);
for t = 0:(N-1)
    idx_up = 3*t + 1;  % pozícia u_t^+
    idx_um = 3*t + 2;  % pozícia u_t^-
    M(t+1, idx_up) =  1;
    M(t+1, idx_um) = -1;
end
% Potom u = M*X a H*(M*X) = x_fin
Aeq = H * M;
beq = x_fin;

% (3c) Nerovnosti pre z_t >= |u(t)| a z_t >= 2|u(t)| - 1
%      kde |u(t)| = u_t^+ + u_t^-:
%
%  (1) z_t >= u_t^+ + u_t^- 
%      ekvivalentne  z_t - (u_t^+ + u_t^-) >= 0
%      čiže          -z_t + u_t^+ + u_t^- <= 0
%
%  (2) z_t >= 2(u_t^+ + u_t^- ) - 1
%      ekvivalentne  z_t - 2(u_t^+ + u_t^-) >= -1
%      čiže          -z_t + 2u_t^+ + 2u_t^- <= 1
%
%  Tieto podmienky zapisujeme v tvare A_ineq * X <= b_ineq.

A_ineq = zeros(2*N, 3*N);
b_ineq = zeros(2*N, 1);

for t = 0:(N-1)
    idx_up = 3*t + 1;  % u_t^+
    idx_um = 3*t + 2;  % u_t^-
    idx_z  = 3*t + 3;  % z_t

    row1 = 2*t + 1;    % prvá nerovnosť
    row2 = 2*t + 2;    % druhá nerovnosť

    % (1) -z_t + u_t^+ + u_t^- <= 0
    A_ineq(row1, idx_z)  = -1; 
    A_ineq(row1, idx_up) =  1; 
    A_ineq(row1, idx_um) =  1; 
    b_ineq(row1)         =  0;

    % (2) -z_t + 2*u_t^+ + 2*u_t^- <= 1
    A_ineq(row2, idx_z)  = -1;
    A_ineq(row2, idx_up) =  2;
    A_ineq(row2, idx_um) =  2;
    b_ineq(row2)         =  1;
end

%% (3d) Dolné a horné hranice pre X
lb = zeros(3*N,1);  % všetky premenné >= 0
ub = inf(3*N,1);    % bez hornej hranice

%% (4) RIEŠENIE POMOCOU linprog
options = optimoptions('linprog','Display','iter','Algorithm','dual-simplex');

[X_opt, fval, exitflag, output] = linprog(c, A_ineq, b_ineq, Aeq, beq, lb, ub, options);

fprintf('----------------------------------------------------\n');
fprintf(' Stav riesenia (exitflag) = %d\n', exitflag);
fprintf(' Minimalna hodnota cielovej funkcie F = %g\n', fval);
fprintf(' Pocet iteracii = %d\n', output.iterations);

%% (5) OBNOVENIE u(t) A z(t)
u_opt = zeros(N,1);
z_opt = zeros(N,1);
for t = 0:(N-1)
    idx_up = 3*t + 1;
    idx_um = 3*t + 2;
    idx_z  = 3*t + 3;
    u_opt(t+1) = X_opt(idx_up) - X_opt(idx_um);
    z_opt(t+1) = X_opt(idx_z);
end

%% (6) (Volitelne) VYPOCET x(t) PRE KAZDY t
x_vals = zeros(n, N+1);
x_vals(:,1) = zeros(n,1);  % x(0)=0
for t = 1:N
    x_vals(:,t+1) = A_mat * x_vals(:,t) + b_vec * u_opt(t);
end

%% (7) VYPIS VYSLEDKOV A ICH GRAFICKE ZOBRAZENIE
disp('Optimalne riadenie u(t):');
disp(u_opt);

disp('Overenie x(N) =');
disp(x_vals(:,end));

figure;
subplot(2,1,1);
stem(0:N-1, u_opt, 'LineWidth',1.5);
title('Optimálne u(t)');
xlabel('t');
ylabel('u(t)');

subplot(2,1,2);
plot(0:N, x_vals','o-','LineWidth',1.3);
legend('x_1(t)','x_2(t)','x_3(t)','Location','Best');
title('Priebeh stavu x(t)');
xlabel('t');
ylabel('x_i(t)');
grid on;
