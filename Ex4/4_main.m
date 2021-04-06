%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Exercice 4  Kefan S., Hiba S., Victor K. Vinh N.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
clc;

%% Scenario 1

disp('===   Scenario 1   ===')

% Création du jeu de données
N = 15; % car il y a : - 1 centrale à charbon
        %              - 5 éoliennes
        %              - 1 barrage
        %              - 2 panneaux photovoltaïques
        %              - 1 data center
        %              - 1 ensemble de 7500 logements
        %              - 1 usine
        %              - 2 trameway
        %              - 1 hôpital
p0  = [ 8.5; 1.9*ones(5,1);   3; 1.9*ones(2,1); -10; -7.5; -9; -.12*ones(2,1); -.2];
a   = [  10;     ones(5,1); .1 ;  .9*ones(2,1);   5;    1;  5;    6*ones(2,1); 200];
b   = [1000;  10*ones(5,1); 100;  11*ones(2,1);  20;   10;  8;    9*ones(2,1);  30];
P_max = [10;   2*ones(5,1);  16;   2*ones(2,1); zeros(6,1)];

% Écriture du problèmean
A   = sparse(1:N,1:N,2*a,N,N);
B   = 2*a.*p0;
C_e = ones(1,N);
d_e = 0;
C_i = sparse([eye(N);[-eye(9) zeros(9,6)]]);
d_i = [P_max ; zeros(9,1)];
cste = sum(a.* p0.^2 + b);

%== Decomposition par prix
lam0 = 0;
mu0  = zeros(N+9,1);
rho  = 0.1;
tol  = 10^-4;
[ U_p, Lam_p, Mu_p, k_p ] = decomp_prix(A,B,C_e,d_e,C_i,d_i, lam0,mu0, rho,tol);

format shortG % format d'affichage
disp('# Décomposition par prix :')
disp(['  - le prix est ', num2str(Lam_p(k_p))]);
cost_p = 1/2* U_p(:,k_p)'*A*U_p(:,k_p) - B'*U_p(:,k_p) + cste;
disp(['  - le coût est ', num2str(cost_p)]);
disp('  - et les puissances fournies sont :')
disp(U_p(:,k_p)')
disp('')


%== Decomposition par ressources
v0  = [d_e;d_i]/N*ones(1,N);
rho = .03;
tol = 10^-4;
rho_uzawa = .1;
tol_uzawa = 10^-4;
[ u_r, v_r, lam_r, mu_r, k_r ] = decomp_ressources_ex4( A,B,C_e,d_e,C_i,d_i, v0, rho,tol, rho_uzawa,tol_uzawa );

disp('# Décomposition par ressources :')
disp(['  - le prix est : ', num2str(mean(lam_r(1,:)))]);
cost_r = 1/2* u_r'*A*u_r - B'*u_r + cste;
disp(['  - le coût est : ', num2str(cost_r)]);
disp('  - et les puissances fournies sont :')
disp(u_r')

%--------------------------------------------------------------------------------
%% Scenario 2
% ATTENTION : ce bout de code met plus d'une minute à s'exécuter

disp('')
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
disp('=== Scenario 2 ===')

% Modification du jeu de données : la puissance maximale des éoliennes est maintenant nulle
P_max = [10; zeros(5,1); 16; 2*ones(2,1); zeros(6,1)];
d_i = [P_max ; zeros(9,1)];


%== Decomposition par prix
lam0 = 0;
mu0  = zeros(N+9,1);
rho  = 0.1;
tol  = 10^-4;
[ U_p, Lam_p, Mu_p, k_p ] = decomp_prix(A,B,C_e,d_e,C_i,d_i, lam0,mu0, rho,tol);

format shortG % format d'affichage
disp('# Décomposition par prix :')
disp(['  - le prix est ', num2str(Lam_p(k_p))]);
cost_p = 1/2* U_p(:,k_p)'*A*U_p(:,k_p) - B'*U_p(:,k_p) + cste;
disp(['  - le coût est ', num2str(cost_p)]);
disp('  - et les puissances fournies sont :')
disp(U_p(:,k_p)')
disp('')


%== Decomposition par ressources
v0  = [d_e;d_i]/N*ones(1,N);
rho = .005;
tol = 10^-3;
rho_uzawa = .02;
tol_uzawa = 10^-3;
[ u_r, v_r, lam_r, mu_r, k_r ] = decomp_ressources_ex4( A,B,C_e,d_e,C_i,d_i, v0, rho,tol, rho_uzawa,tol_uzawa );

disp('# Décomposition par ressources :')
disp(['  - le prix est : ', num2str(mean(lam_r(1,:)))]);
cost_r = 1/2* u_r'*A*u_r - B'*u_r + cste;
disp(['  - le coût est : ', num2str(cost_r)]);
disp('  - et les puissances fournies sont :')
disp(u_r')
