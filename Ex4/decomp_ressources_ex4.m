function [ u_k, v_k, lam_k, mu_k, k ] = decomp_ressources_ex4( A,b,C_e,d_e,C_i,d_i, v0, rho,tol, rho_uzawa,tol_uzawa )
% Algorthime de d�compostion par les ressources
% pour un probl�me quadratique de type :
%   minimiser  J(u) := 1/2<A*u,u> - <b,u>
%      sous     C_e*u == d_e        (e pour contraintes d'�galit�)
%               C_i*u <= d_i        (i pour contraintes d'in�galit�)
%
% Les autres param�tre sont :
% - v0 : une allocation de ressources initiale
% - rho : pas de convergence � l'�tape de coordination (mise � jour des ressources)
% - tol : tol�rance pour l'erreur de convergence
% - rho_uzawa,tol_uzawa : pas et tol�rance pour l'algorithme d'Uzawa utilis�
%
% La fonction renvoie :
% - la solution u_k du probl�me obtenue, 
% - les allocations de ressources v_k associ�es,
% - les multiplicateurs de Lagrange associ�s lam_k et/ou mu_k,
% - ainsi que le nombre d'it�rations effectu�

iter_max = 5000;
N = length(b); % nombre de sous-probl�mes

u_k = zeros(N,1);
v_k = v0; % de dimension m x N, en colonne l'allocation de chaque sous-probl�me (m nombre de contraintes)     (par exemple d/N*ones(1,N))

lam_k = zeros(length(d_e),N);
mu_k = zeros(length(d_i),N);

m = length(d_e) + length(d_i); % nombre de contraintes 
nb_ega = length(d_e); % nombre de contraintes d'�galite, il y a donc m - nb_ega contraintes d'in�galit�
v_k_e = v_k(1:nb_ega,:); % s�paration des allocations li�es aux contraintes d'�galit� et d'in�galit�
v_k_i = v_k((nb_ega+1):m,:);

% Premi�re it�ration (k = 1)

%  i) D�composition :
for i = 1:N
    [u_k(i),lam_k(:,i),mu_k(:,i),~] = Uzawa( A(i,i),b(i),C_e(:,i),v_k_e(:,i),C_i(:,i),v_k_i(:,i), lam_k(:,i),mu_k(:,i), rho_uzawa,tol_uzawa);
end
%  ii) Coordination
v_k_e = v_k_e + rho*(lam_k - 1/N*sum(lam_k,2)*ones(1,N)); % *ones(1,N) permet de faire le calcul de mani�re matriciel
v_k_i = v_k_i + rho*(mu_k - 1/N*sum(mu_k,2)*ones(1,N));

err = tol + 1;

k = 1;
while err > tol && k < iter_max
    k = k+1; 
    u_old = u_k;
    v_old = v_k;
    
    % i) D�composition :
    for i = 1:N
        [u_k(i),lam_k(:,i),mu_k(:,i),~] = Uzawa(A(i,i),b(i),C_e(:,i),v_k_e(:,i),C_i(:,i),v_k_i(:,i), lam_k(:,i),mu_k(:,i), rho_uzawa,tol_uzawa);
    end
    % ii) Coordination
    v_k_e = v_k_e + rho*(lam_k - 1/N*sum(lam_k,2)*ones(1,N)); % *ones(1,N) permet de faire le calcul de mani�re matriciel
    v_k_i = v_k_i + rho*(mu_k - 1/N*sum(mu_k,2)*ones(1,N));
    
    v_k = [v_k_e;v_k_i];
    err = norm(u_k - u_old) + norm(v_k - v_old);
end



if k == iter_max
    disp(['ATTENTION : nombre maximal d''it�rations (',num2str(iter_max),') atteint dans ''decomp_ressources_ex4'''])
end

end