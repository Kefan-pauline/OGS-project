function [ U, Lam, Mu, k ] = DecompPrix(A,b,C_e,d_e,C_i,d_i, lam0,mu0,rho,tol)
% Decomposition par prix pour un problème quadratique de type :
%   minimiser  J(u) := 1/2<A*u,u> - <b,u>
%      sous     C_e*u == d_e        (e pour contraintes d'égalité)
%               C_i*u <= d_i        (i pour contraintes d'inégalité)
%
% La fonction renvoie :
% - la suite U = (u_k) qui doit converger vers la solution du problème, 
% - les multiplicateurs de Lagrange associés Lam = (lam_k) et/ou Mu = (mu_k),
% - ainsi que le nombre d'itérations effectué
iter_max = 1000; % nombre maximal d'itération autorisé
U = zeros(length(b),1); % contiendra au fur à mesure des itérations les valeurs u_k, u_k+1, etc.
Lam = lam0; % intialisation 
Mu = mu0; % intialisation 
N = length(b); % nombre de sous-problèmes

% Première itération
% i) Decomposition
for i=1:N
    U(i,1) = A(i,i)\(b(i) - C_e(:,i)' * Lam(:,1) - C_i(:,i)' * Mu(:,1));
end
% ii) Coordination
Lam(:,2) = Lam(:,1) + rho*(C_e*U(:,1) - d_e);
Mu(:,2) = max(0, Mu(:,1) + rho*(C_i*U(:,1) - d_i));

err = tol + 1; % initialisation de l'erreur (supérieur à la tolérance)
k = 2; % compteur du nombre d'itérations
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
    disp('ATTENTION : nombre maximal d''itérations atteint dans ''Uzawa''')
end
end
