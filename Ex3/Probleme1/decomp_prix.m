function [ U, Lam, Mu, k ] = decomp_prix( A,b,C_e,d_e,C_i,d_i, lam0,mu0, rho,tol )
% Algorthime de décompostion par les prix
% pour un problème quadratique de type (ici l'algorithme est équivalent à Uzawa) :
%   minimiser  J(u) := 1/2<A*u,u> - <b,u>        ( A doit être une matrice symétrique )
%      sous     C_e*u == d_e        (e pour contraintes d'égalité)
%               C_i*u <= d_i        (i pour contraintes d'inégalité)
%
% Les autres paramètre sont :
% - lam0 : un prix initial pour les contraintes d'égalité (mettre 0 si inexistantes)
% - mu0 : un prix initial pour les contraintes d'inégalité (mettre 0 si inexistantes) 
% - rho : pas de convergence à l'étape de coordination (mise à jour des ressources)
% - tol : tolérance pour l'erreur de convergence
%
% La fonction renvoie :
% - la suite U = (u_k) qui doit converger vers la solution du problème, 
% - les multiplicateurs de Lagrange associés Lam = (lam_k) et/ou Mu = (mu_k),
% - ainsi que le nombre d'itérations effectué

iter_max = 2000;
N = length(b); % nombre de sous-problèmes

U = zeros(N,1);
Lam = lam0; % en colonne, les multiplicateurs de Lagrange associés aux N sous-problèmes
Mu = mu0;

% Dans le cas ou seul des contraintes d'égalité ou d'inégalité ont été entrés
if C_e == 0
  C_e = zeros(1,length(b));
end
if C_i == 0
  C_i = zeros(1,length(b));
end

% Première itération (k = 1)
%  i) Décomposition :
%for i = 1:N
%    U(i,1) = A(i,i)\(b(i) - C_e(:,i)'*Lam(:,1) - C_i(:,i)'*Mu(:,1));
%end
U(:,1) = A\(b - C_e'*Lam(:,1) - C_i'*Mu(:,1));

%  ii) Coordination
Lam(:,2) = Lam(:,1) + rho*(C_e*U(:,1) - d_e);
Mu(:,2) = max(0, Mu(:,1) + rho*(C_i*U(:,1) - d_i));

err = tol + 1;

k = 1;
while err > tol && k < iter_max
    k = k+1; 
    % i) Décomposition :
    %for i = 1:N
    %    U(i,k) = A(i,i)\(b(i) - C_e(:,i)'*Lam(:,k) - C_i(:,i)'*Mu(:,k));
    %end
    U(:,k) = A\(b - C_e'*Lam(:,k) - C_i'*Mu(:,k));
    
    % ii) Coordination
    Lam(:,k+1) = Lam(:,k) + rho*(C_e*U(:,k) - d_e);
    Mu(:,k+1) = max(0, Mu(:,k) + rho*(C_i*U(:,k) - d_i));
    
    err = norm(U(:,k) - U(:,k-1)) + norm(Lam(:,k+1) - Lam(:,k)) + norm(Mu(:,k+1) - Mu(:,k));
end

if k == iter_max
    disp('!!! ATTENTION !!! nombre maximal d''itérations atteint dans ''decomp_prix''')
end

end
