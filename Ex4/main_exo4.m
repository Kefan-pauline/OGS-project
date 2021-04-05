%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Exercice 4  Kefan S., Hiba S., Victor K. Vinh N.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;
% Scenario 1
% Creation du jeu de donnees
N = 15;
a = [10; repelem(1,5)'; 0.1; repelem(0.9,2)'; 5; 1; 5; repelem(6,2)'; 200];
A = sparse(1:N,1:N,2*a,N,N);
p0 = [8.5; repelem(1.9,5)'; 3; repelem(1.9,2)'; -10; -7.5; -9; repelem(-0.12,2)'; -0.2];
b = 2*a.*p0;
b_tilde = [1000; repelem(10,5)'; 100; repelem(11,2)'; 20; 10; 8; repelem(9,2)'; 30];
C_e = repelem (1, N);
d_e = 0;
cst = sum(a.* p0.^2 + b_tilde);
P_max =  [10; repelem(2,5)'; 16; repelem(2,2)'; repelem(0,6)'];
C_i = sparse([eye(15);[-eye(9) zeros(9,6)]]);
d_i = [P_max;repelem(0,9)'];

%%
% Decomposition par prix
rho = 0.1;
lam0 = 0;
mu0= zeros(24,1);
tol=10^(-4);
[ U_p, Lam_p, Mu_p, k_p ] = DecompPrix(A,b,C_e,d_e,C_i,d_i, lam0,mu0,rho,tol);
disp(['Le prix est ',num2str(Lam_p(k_p-1))]);
cost_p=(A*U_p(:,k_p-1))'*U_p(:,k_p-1)/2-b'*U_p(:,k_p-1)+cst;
disp(['Le cout est ',num2str(cost_p)]);
%%
% Decomposition par ressource
v0 = [d_e;d_i]/N*ones(1,N);
rho = 0.03;
rho_uzawa=0.1;
tol=10^(-4);
tol_uzawa = 0.1*tol;
[ U_r, V_r, Lam_r, Mu_r, k_r ] = decomp_ressources( A,b,C_e,d_e,C_i,d_i, v0, rho,tol, rho_uzawa,tol_uzawa );
disp(['Le prix est ',num2str(Lam_r(1,15,k_r-1))]);
cost_r=(A*U_r(:,k_r-1))'*U_r(:,k_r-1)/2-b'*U_r(:,k_r-1)+cst;
disp(['Le cout est ',num2str(cost_r)]);
%%
% Scenario 2
% Creation du jeu de donnees
N = 15;
a = [10; repelem(1,5)'; 0.1; repelem(0.9,2)'; 5; 1; 5; repelem(6,2)'; 200];
A = sparse(1:N,1:N,2*a,N,N);
p0 = [8.5; repelem(1.9,5)'; 3; repelem(1.9,2)'; -10; -7.5; -9; repelem(-0.12,2)'; -0.2];
b = 2*a.*p0;
b_tilde = [1000; repelem(10,5)'; 100; repelem(11,2)'; 20; 10; 8; repelem(9,2)'; 30];
C_e = repelem (1, N);
d_e = 0;
cst = sum(a.* p0.^2 + b_tilde);
P_max =  [10; repelem(0,5)'; 16; repelem(2,2)'; repelem(0,6)'];
C_i = sparse([eye(15);[-eye(9) zeros(9,6)]]);
d_i = [P_max;repelem(0,9)'];
%%
% Decomposition par prix
rho = 0.1;
lam0 = 0;
mu0= zeros(24,1);
tol=10^(-4);
[ U_p, Lam_p, Mu_p, k_p ] = DecompPrix(A,b,C_e,d_e,C_i,d_i, lam0,mu0,rho,tol);
disp(['Le prix est ',num2str(Lam_p(k_p-1))]);
cost_p=(A*U_p(:,k_p-1))'*U_p(:,k_p-1)/2-b'*U_p(:,k_p-1)+cst;
disp(['Le cout est ',num2str(cost_p)]);
%%
% Decomposition par ressource
% Ce bout de code met 28s a executer
v0 = [d_e;d_i]/N*ones(1,N);
rho = 0.002;
rho_uzawa=0.05;
tol=10^(-4);
tol_uzawa = 0.1*tol;
[ U_r, V_r, Lam_r, Mu_r, k_r ] = decomp_ressources( A,b,C_e,d_e,C_i,d_i, v0, rho,tol, rho_uzawa,tol_uzawa );
disp(['Le prix est ',num2str(Lam_r(1,15,k_r-1))]);
cost_r=(A*U_r(:,k_r-1))'*U_r(:,k_r-1)/2-b'*U_r(:,k_r-1)+cst;
disp(['Le cout est ',num2str(cost_r)]);
