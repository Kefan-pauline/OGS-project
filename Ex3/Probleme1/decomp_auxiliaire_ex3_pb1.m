function [ U, J, Lam, Mu, k ] = decomp_auxiliaire_ex3_pb1(e,Q,Re, u0, lam,mu, tol, rho_prix,tol_prix, seuil_KKT )
% Algorthime de r�solution par d�compositon par les prix du probl�me auxiliaire associ� au probl�me 1 de l'exercice 3 
% en consid�rant qu'IL Y A des corr�lations entre les actions :
%    minimiser  1/2 <u,Qu>          o� Q n'est pas (forc�ment) diagonale !!!
%      sous     sum u_i == 1 
%               -<e,u> <= -Re
%                  -u_i <= 0  pour tout i dans 1,...,N 
%
% Les autres param�tres sont :
% - u0 : une valeur intiale pour le vecteur de solution U
% - lam : un prix initial pour la contrainte d'�galit�
% - mu : un prix initial pour les contraintes d'in�galit� (mu_0 et le vecteur mu) 
% - rho : pas de convergence � l'�tape de coordination (mise � jour des ressources)
% - tol : tol�rance pour l'erreur de convergence
% - rho_prix : pas dans la d�composition par prix 
% - tol_prix : tol�rance dans la d�composition par prix 
% - seuil_KKT : seuil de respect des conditions KKT � la convergence
%
% La fonction renvoie :
% - la suite U = (u_k) qui doit converger vers la solution du probl�me, 
% - la valeur optimale trouv�e J,
% - les multiplicateurs de Lagrange associ�s Lam = (lam_k) et Mu = (mu_k),
% - ainsi que le nombre d'it�rations effectu�


iter_max = 5000;
N = length(e); % nombre de sous-probl�mes

U = u0;
Lam = lam; % multiplicateurs de Lagrange associ�s aux contraintes
Mu = mu; 

% �criture du probl�me auxilaire pour qu'il puisse �tre entr� dans l'algorithme de d�composition par le prix
D = diag(diag(Q));
%b = -(Q-2*D)*u_k;    b devra �tre r�d�fini dans chaque boucle
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
    disp('!!! ATTENTION !!! nombre maximal d''it�rations atteint dans ''decomp_auxiliaire_ex3_pb1''')
end

conditions_KKT(U(:,k),Lam(:,k),Mu(:,k), 2*D,b,C_e,d_e,C_i,d_i, seuil_KKT, true);

end
