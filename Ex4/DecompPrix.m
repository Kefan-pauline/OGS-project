function [ U, Lam, Mu, k ] = DecompPrix(A,b,C_e,d_e,C_i,d_i, lam0,mu0,rho,tol)
% Decomposition par prix pour un probl�me quadratique de type :
%   minimiser  J(u) := 1/2<A*u,u> - <b,u>
%      sous     C_e*u == d_e        (e pour contraintes d'�galit�)
%               C_i*u <= d_i        (i pour contraintes d'in�galit�)
%
% La fonction renvoie :
% - la suite U = (u_k) qui doit converger vers la solution du probl�me, 
% - les multiplicateurs de Lagrange associ�s Lam = (lam_k) et/ou Mu = (mu_k),
% - ainsi que le nombre d'it�rations effectu�
iter_max = 1000; % nombre maximal d'it�ration autoris�
U = zeros(length(b),1); % contiendra au fur � mesure des it�rations les valeurs u_k, u_k+1, etc.
Lam = lam0; % intialisation 
Mu = mu0; % intialisation 
N = length(b); % nombre de sous-probl�mes

% Premi�re it�ration
% i) Decomposition
for i=1:N
    U(i,1) = A(i,i)\(b(i) - C_e(:,i)' * Lam(:,1) - C_i(:,i)' * Mu(:,1));
end
% ii) Coordination
Lam(:,2) = Lam(:,1) + rho*(C_e*U(:,1) - d_e);
Mu(:,2) = max(0, Mu(:,1) + rho*(C_i*U(:,1) - d_i));

err = tol + 1; % initialisation de l'erreur (sup�rieur � la tol�rance)
k = 2; % compteur du nombre d'it�rations
while err > tol && k <= iter_max
    % i) decomposition
    for i=1:N
        U(i,k) = A(i,i)\(b(i) - C_e(:,i)' * Lam(:,k) - C_i(:,i)' * Mu(:,k));
    end
    % ii) coordination
    Lam(:,k+1) = Lam(:,k) + rho*(C_e*U(:,k) - d_e);
    Mu(:,k+1) = max(0, Mu(:,k) + rho*(C_i*U(:,k) - d_i));
    err = norm(U(:,k) - U(:,k-1));
    k = k+1;
end
if k == iter_max+1
    disp('ATTENTION : nombre maximal d''it�rations atteint dans ''Uzawa''')
end
end
