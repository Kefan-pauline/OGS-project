%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Exercice 3.1 (Kefan S., Hiba S., Victor K. Vinh N.  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------%
%   Problème 1 : imposer un rendement minimum   %
%-----------------------------------------------%

clear
clc

% Génération de données : on simule des actions par des lois normales
N = 4;

rand('seed', 4);  % Graines pour la génération de nombres aléatoires
randn('seed', 4);
sd = rand(1,N)*3 + 1; % écarts-types entre 1 et 4
m  = rand(1,N)*3 + 10; % moyennes entre 10 et 13
Xi = sd.*randn(12,N) + m;

% Affichage des actions simulées
figure(1)
plot(Xi)
axis([1,12])
set(gca,'xtick',1:12,'xticklabel',{'Jan','Fév','Mar','Avr','Mai','Juin','Juil','Août','Sep','Oct','Nov','Déc'})
title('Actions que l''on étudie')


%%-----------------------------------------------------------------------------------------------------------------
% Cas SANS corrélation entre les actions

disp('')
disp('=====  Exercice 3 : Problème 1 (imposer un rendement minimum SANS corrélation)  =====')
disp('')

e  = mean(Xi)';            % Gain moyen de chaque action
Q  = diag(diag(cov(Xi)));  % Risque encouru pour chaque action
Re = mean(e);              % Rendement minimal attendu

% Écriture du problème quadratique :
A = Q;
b = zeros(N,1);
C_e = ones(1,N);
d_e = 1;
C_i = [-e' ; -eye(N)];
d_i = [-Re ; zeros(N,1)];


%== Résolution avec la décomposition par prix (équivalent à Uzawa dans ce cas quadratique) :
lam0 = 0;
mu0  = zeros(N+1,1);
rho  = .0002;
tol  = 10^-5;

[ U, Lam, Mu, k ] = decomp_prix(A,b,C_e,d_e,C_i,d_i, lam0,mu0, rho,tol);
U_decomp_prix = U(:,k)
valeur_optimale = 1/2* U(:,k)'*Q*U(:,k)
nombre_iterations = k
conditions_KKT(U(:,k),Lam(:,k),Mu(:,k), A,b,C_e,d_e,C_i,d_i, .05, true);

disp('')
disp('- - - - - - - - - - - - - - - - - - - - - - - - -')
disp('')

%== Résolution avec la décomposition par quantités :
v0  = [d_e;d_i]/N*ones(1,N);
rho = 3;
tol = 10^-3;
rho_uzawa = .002;
tol_uzawa = 10^-5;

[ U, V, Lam, Mu, k ] = decomp_ressources( A,b,C_e,d_e,C_i,d_i, v0, rho,tol, rho_uzawa,tol_uzawa );
U_decomp_ressources = U(:,k)
valeur_optimale = 1/2* U(:,k)'*Q*U(:,k)
nombre_iterations = k
conditions_KKT(U(:,k),Lam(:,1,k),Mu(:,1,k), A,b,C_e,d_e,C_i,d_i, .05, true);

disp('')
disp('- - - - - - - - - - - - - - - - - - - - - - - - -')
disp('')

%== Résolution avec la décomposition par prédiction :
i0    = 1;
v0    = zeros(N+2,1);
lam0  = 0;
mu0   = zeros(N+1,1);
beta  = .2;
gamma = 1;
tol   = 10^-5;
rho_uzawa = .0002;
tol_uzawa = 10^-5;

[ U, V, Lam, Mu, k ] = decomp_pred_para( A,b,C_e,d_e,C_i,d_i, i0, v0,lam0,mu0, beta,gamma, tol, rho_uzawa,tol_uzawa );
U_decomp_pred = U(:,k)
valeur_optimale = 1/2* U(:,k)'*Q*U(:,k)
nombre_iterations = k
conditions_KKT(U(:,k),Lam(:,k),Mu(:,k), A,b,C_e,d_e,C_i,d_i, .05, true);


%%-----------------------------------------------------------------------------------------------------------------
% Cas AVEC corrélations entre les actions
% avec uniquement des corrélations positives (valeur absolue de la matrice de covariance) 

disp('')
disp('')
disp('=====  Exercice 3 : Problème 1 (imposer un rendement minimum AVEC corrélations)  =====')
disp(' ---   avec des covariances négatives dans la matrice de covariances   ---')
disp('')

e  = mean(Xi)';    % Gain moyen de chaque action
Q  = cov(Xi);      % Risque encouru pour chaque action
Re = mean(e);      % Rendement minimal autorisé

%== Résolution avec la décomposition par prix sur le problème auxiliaire :
u0 = zeros(N,1);
lam = 0;
mu = zeros(N+1,1);
tol = 10^-4;
rho_prix = .002;
tol_prix = 10^-4;
seuil_KKT = 0.05;

[ U, J, Lam, Mu, k ] = decomp_auxiliaire_ex3_pb1(e,Q,Re, u0, lam,mu, tol, rho_prix,tol_prix, seuil_KKT);
U_auxiliaire_corr_negatives = U(:,k)
valeur_optimale = J
nombre_iterations = k


%%-----------------------------------------------------------------------------------------------------------------
% Cas AVEC corrélations entre les actions
% avec uniquement des corrélations positives (valeur absolue de la matrice de covariance)

disp('')
disp('')
disp('=====  Exercice 3 : Problème 1 (imposer un rendement minimum AVEC corrélations)  =====')
disp(' ---   avec des covariances positives uniquement   ---')
disp('')

e  = mean(Xi)';    % Gain moyen de chaque action
Q  = abs(cov(Xi)); % Risque encouru pour chaque action
Re = mean(e);      % Rendement minimal autorisé

%== Résolution avec la décomposition par prix sur le problème auxiliaire :
u0 = zeros(N,1);
lam = 0;
mu = zeros(N+1,1);
tol = 10^-4;
rho_prix = .0019;
tol_prix = 10^-4;
seuil_KKT = 0.05;

[ U, J, Lam, Mu, k ] = decomp_auxiliaire_ex3_pb1(e,Q,Re, u0, lam,mu, tol, rho_prix,tol_prix, seuil_KKT);
format shortG % Format d'affichage (pas en écriture scientifique)
U_auxiliaire_corr_positives = round(U(:,k)*10000)/10000
valeur_optimale = J
nombre_iterations = k
