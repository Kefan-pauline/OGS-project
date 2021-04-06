function [ U, J, Lam, Mu, k ] = decomp_reformulation_ex3_pb2( e,Q,De,alpha, u0,lam,mu, tol, rho_prix,tol_prix, seuil_KKT)
% Algorthime de d�compostion par les prix de la reformulation du probl�me 2 de l'exercice 3 
% en consid�rant qu'IL Y A des corr�lations entre les actions :
%    maximiser  <e,u>
%      sous     sum u_i == 1 
%                <u,Qu> <= De        o� Q n'est pas (forc�ment) diagonale !!!
%                  -u_i <= 0  pour tout i dans 1,...,N 
%
% reformul� de la mani�re suivante :
%    minimiser  alpha*<u,u> - <e+2*alpha*u_k,u>
%      sous     sum u_i == 1
%             <u,Q*u_k> <= De
%                  -u_i <= 0  pour tout i dans 1,...,N 
%
% Les autres param�tres sont :
% - u0 : une valeur intiale pour le vecteur de solution U
% - lam : un prix initial pour la contrainte d'�galit�
% - mu : un prix initial pour les contraintes d'in�galit� (mu_0 et le vecteur mu) 
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


iter_max = 1000;
N = length(e); % nombre de sous-probl�mes

U = u0;
Lam = lam; % multiplicateurs de Lagrange associ�s aux contraintes
Mu = mu; 

% �criture du probl�me auxilaire pour qu'il puisse �tre entr� dans l'algorithme de d�composition par le prix
A = 2*alpha*eye(N);
%b = e + 2*alpha*u_k;        % b devra �tre r�d�fini dans chaque boucle car d�pend de u_k
C_e = ones(1,N);
d_e = 1;
%C_i = [u_k'*Q ; -eye(N)];   % C_i devra aussi �tre r�d�fini dans chaque boucle car il d�pend aussi de u_k
d_i = [De ; zeros(N,1)];

err = tol + 1;

k = 1;
while err > tol && k < iter_max
    % Mise � jour de b et de C_i :
    b = e + 2*alpha*U(:,k);
    C_i = [U(:,k)'*Q ; -eye(N)];
    
    [ U_dp, Lam_dp, Mu_dp, k_dp ] = decomp_prix( A,b,C_e,d_e,C_i,d_i, Lam(k),Mu(:,k), rho_prix,tol_prix );
    U(:,k+1) = U_dp(:,k_dp);
    Lam(k+1) = Lam_dp(:,k_dp);
    Mu(:,k+1) = Mu_dp(:,k_dp);
    
    err = norm(U(:,k+1) - U(:,k)) + norm([Lam(k+1) - Lam(k); Mu(:,k+1) - Mu(:,k)]);
    k = k+1;
end

J = e'*U(:,k);

if k == iter_max
    disp('!!! ATTENTION !!! nombre maximal d''it�rations atteint dans ''decomp_reformulation_ex3_pb2''')
end

conditions_KKT(U(:,k),Lam(:,k),Mu(:,k), A,b,C_e,d_e,C_i,d_i, seuil_KKT, true);

end
