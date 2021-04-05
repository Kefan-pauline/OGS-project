function [ U, J, Lam, Mu, k ] = decomp_prix_ex3_pb2_sans_corr( e,Q,De, lam,mu, rho,tol, seuil_KKT )
% Algorthime de décompostion par les prix pour le problème 2 de l'exercice 3 
% en considérant qu'IL N'Y A PAS de corrélations entre les actions :
%    maximiser  <e,u>
%      sous     sum u_i == 1 
%                <u,Qu> <= De                     où Q est diagonale !!!
%                  -u_i <= 0  pour tout i dans 1,...,N 
%
% Les autres paramètres sont :
% - lam : un prix initial pour la contrainte d'égalité
% - mu : un prix initial pour les contraintes d'inégalité (mu_0 et le vecteur mu) 
% - rho : pas de convergence à l'étape de coordination (mise à jour des ressources)
% - tol : tolérance pour l'erreur de convergence
% - seuil_KKT : seuil de respect des conditions KKT à la convergence
%
% La fonction renvoie :
% - la suite U = (u_k) qui doit converger vers la solution du problème, 
% - la valeur optimale trouvée J,
% - les multiplicateurs de Lagrange associés Lam = (lam_k) et Mu = (mu_k),
% - ainsi que le nombre d'itérations effectué


iter_max = 10000;
N = length(e); % nombre de sous-problèmes

U = zeros(N,1);
Lam = lam; % multiplicateurs de Lagrange associés aux contraintes
Mu = mu; 

% Écriture des contraintes d'égalité (en haut) et d'inégalité (en-dessous)
Theta_egalite = @(u) sum(u);
Theta_inegalites = @(u) [u'*Q*u ; -u];
d_e = 1;
d_i = [De ; zeros(length(e),1)];

% Première itération (k = 1)
%  i) Décomposition :
for i = 1:N
    U(i,1) = (2*Mu(1,1)*Q(i,i))\(Mu(1+i,1) + e(i) - Lam(1));
end

%  ii) Coordination
Lam(2) = Lam(1) + rho*( Theta_egalite(U(:,1)) - d_e );
Mu(:,2) = max(0, Mu(:,1) + rho*( Theta_inegalites(U(:,1)) - d_i ) );

err = tol + 1;

k = 1;
while err > tol && k < iter_max
    k = k+1; 
    % i) Décomposition :
    for i = 1:N
        U(i,k) = (2*Mu(1,k)*Q(i,i))\(Mu(1+i,k) + e(i) - Lam(k));
    end
    
    % ii) Coordination
    Lam(k+1) = Lam(k) + rho*( Theta_egalite(U(:,k)) - d_e );
    Mu(:,k+1) = max(0, Mu(:,k) + rho*( Theta_inegalites(U(:,k)) - d_i) );
    
    err = norm([Lam(k+1) - Lam(k); Mu(:,k+1) - Mu(:,k)]);
end

J = e'*U(:,k);

if k == iter_max
    disp('!!! ATTENTION !!! nombre maximal d''itérations atteint dans ''decomp_prix_ex3_pb2_sans_corr''')
end

% Conditions KKT :
if all(abs(- e + Lam(k) + 2*Mu(1,k)*Q*U(:,k) - Mu(2:(length(e)+1),k)) < seuil_KKT) ...
   && all(Theta_egalite(U(:,k)) - d_e < seuil_KKT) ...                                  % pour les contraintes d'égalité
   && all(Theta_inegalites(U(:,k)) - d_i < seuil_KKT) && all(Mu(:,k) > -seuil_KKT) ...  % pour les contraintes d'inégalité
   && all(abs(Mu(:,k) .* (Theta_inegalites(U(:,k)) - d_i)) < seuil_KKT); 
     disp(['Les conditions KKT sont vérifiées au seuil d''erreur de ',num2str(seuil_KKT)])
else
     disp(['!!! ATTENTION !!! Les conditions KKT NE sont PAS vérifiées au seuil d''erreur de ',num2str(seuil_KKT)])
end

end
