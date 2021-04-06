%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Exercice 3.1 (Kefan S., Hiba S., Victor K. Vinh N.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------------------------------%
%   Probl�me 1 : imposer un rendement minimum   %
%-----------------------------------------------%

clear
clc

% G�n�ration de donn�es : on simule des actions par des lois normales
N = 4;

%rand('seed', 4);  % Graines pour la g�n�ration de nombres al�atoires
%randn('seed', 4);
%sd = rand(1,N)*3 + 1; % �carts-types entre 1 et 4
%m  = rand(1,N)*3 + 10; % moyennes entre 10 et 13
%Xi = (ones(12,1)*sd).*randn(12,N) + ones(12,1)*m;

% Afin de s'assurer que les valeurs utilis�es sont les m�mes que celle pr�sent�es dans le rapport
% on entre "� la main" les donn�es utilis�es :
sd = [1.3855, 2.7536, 1.8116, 1.5314];
m = [10.709, 12.386, 12.584, 11.748];
Xi =[[11.106, 16.721, 10.618, 12.983],
     [10.468, 12.327, 10.969, 10.716],
     [11.663, 9.1281, 13.305, 9.6623],
     [11.285, 4.6205, 17.734, 13.542],
     [11.517, 14.657, 11.274, 11.53 ],
     [9.6403, 11.699, 17.934, 9.9179],
     [9.3042, 16.305, 9.1041, 10.793],
     [10.475,   10.7, 10.114, 13.505],
     [9.1433, 11.039, 13.506, 11.57 ],
     [  12.8, 12.681, 12.134, 10.689],
     [11.593, 14.207, 11.792, 13.072],
     [10.544, 11.258, 10.518, 10.652]];

     
% Affichage des actions simul�es
figure(1)
plot(Xi)
axis([1,12, 4,19])
set(gca,'xtick',1:12,'xticklabel',{'Jan','F�v','Mar','Avr','Mai','Juin','Juil','Ao�t','Sep','Oct','Nov','D�c'})
title(['Les ',num2str(N),' actions observ�es'])


%%-----------------------------------------------------------------------------------------------------------------
% Cas SANS corr�lation entre les actions

disp('')
disp('=====  Exercice 3 : Probl�me 1 (imposer un rendement minimum SANS corr�lation)  =====')
disp('')

e  = mean(Xi)';            % Gain moyen de chaque action
Q  = diag(diag(cov(Xi)));  % Risque encouru pour chaque action
Re = mean(e);              % Rendement minimal attendu

% �criture du probl�me quadratique :
A = Q;
b = zeros(N,1);
C_e = ones(1,N);
d_e = 1;
C_i = [-e' ; -eye(N)];
d_i = [-Re ; zeros(N,1)];


%== R�solution avec la d�composition par prix (�quivalent � Uzawa dans ce cas quadratique) :
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

%== R�solution avec la d�composition par quantit�s :
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

%== R�solution avec la d�composition par pr�diction :
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
% Cas AVEC corr�lations entre les actions
% avec uniquement des corr�lations positives (valeur absolue de la matrice de covariance) 

disp('')
disp('')
disp('=====  Exercice 3 : Probl�me 1 (imposer un rendement minimum AVEC corr�lations)  =====')
disp(' ---   avec des covariances n�gatives dans la matrice de covariances   ---')
disp('')

e  = mean(Xi)';    % Gain moyen de chaque action
Q  = cov(Xi);      % Risque encouru pour chaque action
Re = mean(e);      % Rendement minimal autoris�

%== R�solution avec la d�composition par prix sur le probl�me auxiliaire :
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
% Cas AVEC corr�lations entre les actions
% avec uniquement des corr�lations positives (valeur absolue de la matrice de covariance)

disp('')
disp('')
disp('=====  Exercice 3 : Probl�me 1 (imposer un rendement minimum AVEC corr�lations)  =====')
disp(' ---   avec des covariances positives uniquement   ---')
disp('')

e  = mean(Xi)';    % Gain moyen de chaque action
Q  = abs(cov(Xi)); % Risque encouru pour chaque action
Re = mean(e);      % Rendement minimal autoris�

%== R�solution avec la d�composition par prix sur le probl�me auxiliaire :
u0 = zeros(N,1);
lam = 0;
mu = zeros(N+1,1);
tol = 10^-4;
rho_prix = .0019;
tol_prix = 10^-4;
seuil_KKT = 0.05;

[ U, J, Lam, Mu, k ] = decomp_auxiliaire_ex3_pb1(e,Q,Re, u0, lam,mu, tol, rho_prix,tol_prix, seuil_KKT);
format shortG % Format d'affichage (pas en �criture scientifique)
U_auxiliaire_corr_positives = round(U(:,k)*10000)/10000
valeur_optimale = J
nombre_iterations = k
