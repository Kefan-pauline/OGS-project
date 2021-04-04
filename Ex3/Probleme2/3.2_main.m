%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Exercice 3.2 (Kefan S., Hiba S., Victor K. Vinh N.  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------%
%   Probl�me 2 : imposer un risque maximum   %
%--------------------------------------------%

clear
clc

% G�n�ration de donn�es : on simule des actions par des lois normales
N = 4;

rand('seed', 4); % Graines pour la g�n�ration de nombres al�atoires
randn('seed', 4);
sd = rand(1,N)*3 + 1; % �carts-types entre 1 et 4
m  = rand(1,N)*3 + 10; % moyennes entre 10 et 13
Xi = sd.*randn(12,N) + m;

% Affichage des actions simul�es
figure(1)
plot(Xi)
axis([1,12, 0,20])
set(gca,'xtick',1:12,'xticklabel',{'Jan','F�v','Mar','Avr','Mai','Juin','Juil','Ao�t','Sep','Oct','Nov','D�c'})
title('Actions que l''on �tudie')


%%-----------------------------------------------------------------------------------------------------------------
% Cas SANS corr�lation entre les actions
% avec des corr�lations n�gatives (vraie matrice de covariance) 

disp('')
disp('=====  Exercice 3 : Probl�me 2 (imposer un risque maximum SANS corr�lation)  =====')
disp('')

e  = mean(Xi)';            % Gain moyen de chaque action
Q  = diag(diag(cov(Xi)));  % Risque encouru pour chaque action
De = min(diag(Q));         % Risque maximal autoris�


%== R�solution avec la d�composition par prix :
lam0 = 0;
mu0 = ones(N+1,1);
rho = .4;
tol = 10^-5;
seuil_KKT = 0.001;

[ U, J, Lam, Mu, k ] = decomp_prix_ex3_pb2_sans_corr( e,Q,De, lam0,mu0, rho,tol, seuil_KKT );
U_decomp_prix_sans_corr = U(:,k)
valeur_optimale = J
nombre_iteration = k



%%-----------------------------------------------------------------------------------------------------------------
% Cas AVEC corr�lations entre les actions
% avec uniquement des corr�lations positives (valeur absolue de la matrice de covariance) 

disp('')
disp('')
disp('=====  Exercice 3 : Probl�me 2 (imposer un risque maximum AVEC corr�lations)  =====')
disp(' ---   avec des covariances n�gatives dans la matrice de covariances   ---')
disp('')

e  = mean(Xi)';    % Gain moyen de chaque action
Q  = cov(Xi);      % Risque encouru pour chaque action (uniquement de covariances postives ici !!!! )
De = min(diag(Q)); % Risque maximal autoris�

% L'�criture du probl�me quadratique se fait directment dans les fonctions car le vecteur b
% ainsi que les contraintes d'in�galit�s doivent �tre remise � jour � chaque it�ration


%== R�solution avec la d�composition par prix de la reformulation du probl�me :
alpha = 3.4;
u0    = zeros(N,1);
lam   = 0;
mu    = zeros(N+1,1);
tol   = 10^-3;
rho_prix = .1;
tol_prix = 10^-4;
seuil_KKT = 0.001;

[ U, J, Lam, Mu, k ] = decomp_reformulation_ex3_pb2( e,Q,De,alpha, u0,lam,mu, tol, rho_prix,tol_prix, seuil_KKT );
format shortG % Format d'affichage (pas en �criture scientifique)
U_decomp_prix_corr_negatives = round(U(:,k)*10000)/10000
valeur_optimale = J
nombre_iteration = k



%%-----------------------------------------------------------------------------------------------------------------
% Cas AVEC corr�lations entre les actions
% avec uniquement des corr�lations positives (valeur absolue de la matrice de covariance) 

disp('')
disp('')
disp('=====  Exercice 3 : Probl�me 2 (imposer un risque maximum AVEC corr�lations)  =====')
disp(' ---   avec des covariances positives uniquement   ---')
disp('')

e  = mean(Xi)';    % Gain moyen de chaque action
Q  = abs(cov(Xi)); % Risque encouru pour chaque action (uniquement de covariances postives ici !!!! )
De = min(diag(Q)); % Risque maximal autoris�

% L'�criture du probl�me quadratique se fait directment dans les fonctions car le vecteur b
% ainsi que les contraintes d'in�galit�s doivent �tre remise � jour � chaque it�ration


%== R�solution avec la d�composition par prix de la reformulation du probl�me :
alpha = 5;
u0    = zeros(N,1);
lam   = 0;
mu    = zeros(N+1,1);
tol   = 10^-3;
rho_prix = .1;
tol_prix = 10^-4;
seuil_KKT = 0.001;

[ U, J, Lam, Mu, k ] = decomp_reformulation_ex3_pb2( e,Q,De,alpha, u0,lam,mu, tol, rho_prix,tol_prix, seuil_KKT );
U_decomp_prix_corr_positives = round(U(:,k)*10000)/10000
valeur_optimale = J
nombre_iteration = k