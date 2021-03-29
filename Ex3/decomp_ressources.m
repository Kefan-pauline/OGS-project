function [ U, V, Lam, Mu, k ] = decomp_ressources( A,b,C_e,d_e,C_i,d_i, v0, rho,tol, rho_uzawa,tol_uzawa )
%% Algorthime de d�compostion par les ressources
%% pour un probl�me quadratique de type :
%%   minimiser  J(u) := 1/2<A*u,u> - <b,u>
%%      sous     C_e*u == d_e        (e pour contraintes d'�galit�)
%%               C_i*u <= d_i        (i pour contraintes d'in�galit�)
%%
%% Les autres param�tre sont :
%% - v0 : une allocation de ressources initiale
%% - rho : pas de convergence � l'�tape de coordination (mise � jour des ressources)
%% - tol : tol�rance pour l'erreur de convergence
%% - rho_uzawa,tol_uzawa : pas et tol�rance pour l'algorithme d'Uzawa utilis�
%%
%% La fonction renvoie :
%% - la suite U = (u_k) qui doit converger vers la solution du probl�me, 
%% - la suite V = (v_k) des allocations de ressources,
%% - les multiplicateurs de Lagrange associ�s Lam = (lam_k) et/ou Mu = (mu_k),
%% - ainsi que le nombre d'it�rations effectu�


iter_max = 1000;
N = length(b); % nombre de sous-probl�mes

U = zeros(N,1);
V = v0; % de dimension m x N, en colonne l'allocation de chaque sous-probl�me (m nombre de contraintes)     (par exemple d/N*ones(1,N))

Lam = zeros(length(d_e),N); % en colonne, les multiplicateurs de Lagrange associ�s aux N sous-probl�mes
Mu = zeros(length(d_i),N);

%-------------------------------------------------------------------------------------------------------------------------------------------
% Cas 1 : s'il n'y a que des contraintes d'�galit�
if C_i == 0
  % Premi�re it�ration (k = 1)
  
  %  i) D�composition :
  for i = 1:N
      [U(i,1),Lam(:,i,1),~,~] = Uzawa( A(i,i),b(i),C_e(:,i),V(:,i,1),0,0, Lam(:,i,1),0, rho_uzawa,tol_uzawa);
  end
  %  ii) Coordination
  V(:,:,2) = V(:,:,1) + rho*(Lam(:,:,1) - 1/N*sum(Lam(:,:,1),2)*ones(1,N)); % *ones(1,N) permet de faire le calcul de mani�re matriciel

  err = tol + 1;

  k = 1;
  while err > tol && k < iter_max
      k = k+1; 
      % i) D�composition :
      for i = 1:N
          [U(i,k),Lam(:,i,k),~,~] = Uzawa(A(i,i),b(i),C_e(:,i),V(:,i,k),0,0, Lam(:,i,k-1),0, rho_uzawa,tol_uzawa);
      end
      % ii) Coordination
      V(:,:,k+1) = V(:,:,k) + rho*(Lam(:,:,k) - 1/N*sum(Lam(:,:,k),2)*ones(1,N));
      
      err = norm(U(:,k) - U(:,k-1));
  end


%-------------------------------------------------------------------------------------------------------------------------------------------
% Cas 2 : s'il n'y a que des contraintes d'in�galit�
elseif C_e == 0 
  % Premi�re it�ration (k = 1)
  
  %  i) D�composition :
  for i = 1:N
      [U(i,1),~,Mu(:,i,1),~] = Uzawa(A(i,i),b(i),0,0,C_i(:,i),V(:,i,1), 0,Mu(:,i,1), rho_uzawa,tol_uzawa);
  end
  %  ii) Coordination :
  V(:,:,2) = V(:,:,1) + rho*(Mu(:,:,1) - 1/N*sum(Mu(:,:,1),2)*ones(1,N)); % *ones(1,N) permet de faire le calcul de mani�re matriciel

  err = tol + 1;

  k = 1;
  while err > tol && k < iter_max
      k = k+1; 
      % i) D�composition :
      for i = 1:N
          [U(i,k),~,Mu(:,i,k),~] = Uzawa(A(i,i),b(i),0,0,C_i(:,i),V(:,i,k), 0,Mu(:,i,k-1), rho_uzawa,tol_uzawa);
      end
      % ii) Coordination
      V(:,:,k+1) = V(:,:,k) + rho*(Mu(:,:,k) - 1/N*sum(Mu(:,:,k),2)*ones(1,N));
      
      err = norm(U(:,k) - U(:,k-1));
  end
  
  
%-------------------------------------------------------------------------------------------------------------------------------------------
% Cas 3 : s'il y a des contraintes d'�galit� et d'in�galit�
else
  m = length(d_e) + length(d_i); # nombre de contraintes 
  nb_ega = length(d_e); # nombre de contraintes d'�galite, il y a donc m - nb_ega contraintes d'in�galit�

  V_e = V(1:nb_ega,:); % s�paration des allocations li�es aux contraintes d'�galit� et d'in�galit�
  V_i = V((nb_ega+1):m,:);
  
  % Premi�re it�ration (k = 1)
  
  %  i) D�composition :
  for i = 1:N
      [U(i,1),Lam(:,i,1),Mu(:,i,1),~] = Uzawa( A(i,i),b(i),C_e(:,i),V_e(:,i,1),C_i(:,i),V_i(:,i,1), Lam(:,i,1),Mu(:,i,1), rho_uzawa,tol_uzawa);
  end
  %  ii) Coordination
  V_e(:,:,2) = V_e(:,:,1) + rho*(Lam(:,:,1) - 1/N*sum(Lam(:,:,1),2)*ones(1,N)); % *ones(1,N) permet de faire le calcul de mani�re matriciel
  V_i(:,:,2) = V_i(:,:,1) + rho*(Mu(:,:,1) - 1/N*sum(Mu(:,:,1),2)*ones(1,N));

  err = tol + 1;

  k = 1;
  while err > tol && k < iter_max
      k = k+1; 
      % i) D�composition :
      for i = 1:N
          [U(i,k),Lam(:,i,k),Mu(:,i,k),~] = Uzawa(A(i,i),b(i),C_e(:,i),V_e(:,i,k),C_i(:,i),V_i(:,i,k), Lam(:,i,k-1),Mu(:,i,k-1), rho_uzawa,tol_uzawa);
      end
      % ii) Coordination
      V_e(:,:,k+1) = V_e(:,:,k) + rho*(Lam(:,:,k) - 1/N*sum(Lam(:,:,k),2)*ones(1,N)); % *ones(1,N) permet de faire le calcul de mani�re matriciel
      V_i(:,:,k+1) = V_i(:,:,k) + rho*(Mu(:,:,k) - 1/N*sum(Mu(:,:,k),2)*ones(1,N));
      
      err = norm(U(:,k) - U(:,k-1));
  end
  V = [V_e;V_i];
  
end


if k == iter_max
    disp('ATTENTION : nombre maximal d''it�rations atteint dans ''decomp_ressources''')
end

end
