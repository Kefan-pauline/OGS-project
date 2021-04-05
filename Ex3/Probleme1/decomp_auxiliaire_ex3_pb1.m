function [ U, J, Lam, Mu, k ] = decomp_auxiliaire_ex3_pb1(e,Q,Re, u0, lam,mu, tol, rho_prix,tol_prix, seuil_KKT )
% Algorthime de résolution par décompositon par les prix du problème auxiliaire associé au problème 1 de l'exercice 3 
% en considérant qu'IL Y A des corrélations entre les actions :
%    minimiser  1/2 <u,Qu>          où Q n'est pas (forcément) diagonale !!!
%      sous     sum u_i == 1 
%               -<e,u> <= -Re
%                  -u_i <= 0  pour tout i dans 1,...,N 
%
% Les autres paramètres sont :
% - u0 : une valeur intiale pour le vecteur de solution U
% - lam : un prix initial pour la contrainte d'égalité
% - mu : un prix initial pour les contraintes d'inégalité (mu_0 et le vecteur mu) 
% - rho : pas de convergence à l'étape de coordination (mise à jour des ressources)
% - tol : tolérance pour l'erreur de convergence
% - rho_prix : pas dans la décomposition par prix 
% - tol_prix : tolérance dans la décomposition par prix 
% - seuil_KKT : seuil de respect des conditions KKT à la convergence
%
% La fonction renvoie :
% - la suite U = (u_k) qui doit converger vers la solution du problème, 
% - la valeur optimale trouvée J,
% - les multiplicateurs de Lagrange associés Lam = (lam_k) et Mu = (mu_k),
% - ainsi que le nombre d'itérations effectué


iter_max = 5000;
N = length(e); % nombre de sous-problèmes

U = u0;
Lam = lam; % multiplicateurs de Lagrange associés aux contraintes
Mu = mu; 

% Écriture du problème auxilaire pour qu'il puisse être entré dans l'algorithme de décomposition par le prix
D = diag(diag(Q));
%b = -(Q-2*D)*u_k;    b devra être rédéfini dans chaque boucle
C_e = ones(1,N);
d_e = 1;
C_i = [-e' ; -eye(N)];
d_i = [-Re ; zeros(N,1)];

err = tol + 1;

k = 1;
while err > tol && k < iter_max
    b = -(Q-2*D)*U(:,k);
    
    [ U_dp, Lam_dp, Mu_dp, k_dp ] = decomp_prix( 2*D,b,C_e,d_e,C_i,d_i, Lam(k),Mu(:,k), rho_prix,tol_prix );
    U(:,k+1) = U_dp(:,k_dp);
    Lam(k+1) = Lam_dp(:,k_dp);
    Mu(:,k+1) = Mu_dp(:,k_dp);
    
    err = norm(U(:,k+1) - U(:,k)) + norm([Lam(k+1) - Lam(k); Mu(:,k+1) - Mu(:,k)]);
    k = k+1; 
end

J = 1/2* U(:,k)'*Q*U(:,k);

if k == iter_max
    disp('!!! ATTENTION !!! nombre maximal d''itérations atteint dans ''decomp_auxiliaire_ex3_pb1''')
end

conditions_KKT(U(:,k),Lam(:,k),Mu(:,k), 2*D,b,C_e,d_e,C_i,d_i, seuil_KKT, true);

end
