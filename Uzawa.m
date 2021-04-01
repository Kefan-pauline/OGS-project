function [ u_k, lam_k, mu_k, k ] = Uzawa( A,b,C_e,d_e,C_i,d_i, lam0,mu0, rho,tol )
%% Algorthime d'Uzawa pour un probl�me quadratique de type :
%%   minimiser  J(u) := 1/2<A*u,u> - <b,u>
%%      sous     C_e*u == d_e        (e pour contraintes d'�galit�)
%%               C_i*u <= d_i        (i pour contraintes d'in�galit�)
%%
%% Les autres param�tres sont :
%% - lam0 : un prix initial pour les contraintes d'�galit� (mettre 0 si inexistantes)
%% - mu0 : un prix initial pour les contraintes d'in�galit� (mettre 0 si inexistantes) 
%% - rho : pas de convergence pour la mise � jour des prix
%% - tol : tol�rance pour l'erreur de convergence
%%
%% La fonction renvoie :
%% - la solution approch�e u_k du probl�me, 
%% - les multiplicateurs de Lagrange associ�s lam_k et mu_k,
%% - ainsi que le nombre d'it�rations effectu�

iter_max = 1000; % nombre maximal d'it�ration autoris�

u_k = zeros(length(b),1); % contiendra au fur � mesure des it�rations les valeurs u_k, u_k+1, etc.
lam_k = lam0; % intialisation avec lam_0 
mu_k = mu0; % intialisation avec mu_0

% Dans le cas ou seul des contraintes d'�galit� ou d'in�galit� ont �t� entr�s
if C_e == 0
  C_e = zeros(1,length(b));
end
if C_i == 0
  C_i = zeros(1,length(b));
end

% Premi�re it�ration (k = 1)
u_k = A\(b - C_e' * lam_k - C_i' * mu_k);
lam_k = lam_k + rho*(C_e*u_k - d_e);
mu_k = max(0, mu_k + rho*(C_i*u_k - d_i));
err = tol + 1; % initialisation de l'erreur (sup�rieur � la tol�rance)

k = 1; % compteur du nombre d'it�rations
while err > tol && k < iter_max
    k = k+1;
    u_old = u_k;
    
    u_k = A\(b - C_e' * lam_k - C_i' * mu_k);
    lam_k = lam_k + rho*(C_e*u_k - d_e);
    mu_k = max(0, mu_k + rho*(C_i*u_k - d_i));
    err = norm(u_k - u_old);
end

if k == iter_max
    disp('ATTENTION : nombre maximal d''it�rations atteint dans ''Uzawa''')
end

end


%{
M�me fonction mais qui garde les U, Lam et Mu en m�moire et les retourne

function [ U, Lam, Mu, k ] = Uzawa(A,b,C_e,d_e,C_i,d_i, rho_e,rho_i,tol)
%% Algorthime d'Uzawa pour un probl�me quadratique de type :
%%   minimiser  J(u) := 1/2<A*u,u> - <b,u>
%%      sous     C_e*u == d_e        (e pour contraintes d'�galit�)
%%               C_i*u <= d_i        (i pour contraintes d'in�galit�)
%%
%% La fonction renvoie :
%% - la suite U = (u_k) qui doit converger vers la solution du probl�me, 
%% - les multiplicateurs de Lagrange associ�s Lam = (lam_k) et/ou Mu = (mu_k),
%% - ainsi que le nombre d'it�rations effectu�

iter_max = 1000; % nombre maximal d'it�ration autoris�

U = zeros(length(b),1); % contiendra au fur � mesure des it�rations les valeurs u_k, u_k+1, etc.
Lam = zeros(length(d_e),1); % intialisation avec lam_0 = (0,...,0)
Mu = zeros(length(d_i),1); % intialisation avec mu_0 = (0,...,0)

% Dans le cas ou seul des contraintes d'�galit� ou d'in�galit� ont �t� entr�s
if C_e == 0
  C_e = zeros(1,length(b));
end
if C_i == 0
  C_i = zeros(1,length(b));
end

% Premi�re it�ration
U(:,1) = A\(b - C_e' * Lam(:,1) - C_i' * Mu(:,1));
Lam(:,2) = Lam(:,1) + rho_e*(C_e*U(:,1) - d_e);
Mu(:,2) = max(0, Mu(:,1) + rho_i*(C_i*U(:,1) - d_i));
err = tol + 1; % initialisation de l'erreur (sup�rieur � la tol�rance)

k = 2; % compteur du nombre d'it�rations
while err > tol && k <= iter_max
    U(:,k) = A\(b - C_e' * Lam(:,k) - C_i' * Mu(:,k));
    Lam(:,k+1) = Lam(:,k) + rho_e*(C_e*U(:,k) - d_e);
    Mu(:,k+1) = max(0, Mu(:,k) + rho_i*(C_i*U(:,k) - d_i));
    err = norm(U(:,k) - U(:,k-1));
    k = k+1;
end

if k == iter_max+1
    disp('ATTENTION : nombre maximal d''it�rations atteint dans ''Uzawa''')
end

end
%}