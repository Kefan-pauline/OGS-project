function [ U, J, Lam, Mu, k ] = decomp_prix_ex3_pb2_sans_corr( e,Q,De, lam,mu, rho,tol, seuil_KKT )
% Algorthime de d�compostion par les prix pour le probl�me 2 de l'exercice 3 
% en consid�rant qu'IL N'Y A PAS de corr�lations entre les actions :
%    maximiser  <e,u>
%      sous     sum u_i == 1 
%                <u,Qu> <= De                     o� Q est diagonale !!!
%                  -u_i <= 0  pour tout i dans 1,...,N 
%
% Les autres param�tres sont :
% - lam : un prix initial pour la contrainte d'�galit�
% - mu : un prix initial pour les contraintes d'in�galit� (mu_0 et le vecteur mu) 
% - rho : pas de convergence � l'�tape de coordination (mise � jour des ressources)
% - tol : tol�rance pour l'erreur de convergence
% - seuil_KKT : seuil de respect des conditions KKT � la convergence
%
% La fonction renvoie :
% - la suite U = (u_k) qui doit converger vers la solution du probl�me, 
% - la valeur optimale trouv�e J,
% - les multiplicateurs de Lagrange associ�s Lam = (lam_k) et Mu = (mu_k),
% - ainsi que le nombre d'it�rations effectu�


iter_max = 10000;
N = length(e); % nombre de sous-probl�mes

U = zeros(N,1);
Lam = lam; % multiplicateurs de Lagrange associ�s aux contraintes
Mu = mu; 

% �criture des contraintes d'�galit� (en haut) et d'in�galit� (en-dessous)
Theta_egalite = @(u) sum(u);
Theta_inegalites = @(u) [u'*Q*u ; -u];
d_e = 1;
d_i = [De ; zeros(length(e),1)];

% Premi�re it�ration (k = 1)
%  i) D�composition :
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
    % i) D�composition :
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
    disp('!!! ATTENTION !!! nombre maximal d''it�rations atteint dans ''decomp_prix_ex3_pb2_sans_corr''')
end

% Conditions KKT :
if all(abs(- e + Lam(k) + 2*Mu(1,k)*Q*U(:,k) - Mu(2:(length(e)+1),k)) < seuil_KKT) ...
   && all(Theta_egalite(U(:,k)) - d_e < seuil_KKT) ...                                  % pour les contraintes d'�galit�
   && all(Theta_inegalites(U(:,k)) - d_i < seuil_KKT) && all(Mu(:,k) > -seuil_KKT) ...  % pour les contraintes d'in�galit�
   && all(abs(Mu(:,k) .* (Theta_inegalites(U(:,k)) - d_i)) < seuil_KKT); 
     disp(['Les conditions KKT sont v�rifi�es au seuil d''erreur de ',num2str(seuil_KKT)])
else
     disp(['!!! ATTENTION !!! Les conditions KKT NE sont PAS v�rifi�es au seuil d''erreur de ',num2str(seuil_KKT)])
end

end
