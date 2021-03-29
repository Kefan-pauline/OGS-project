function [ U, V, Lam, Mu, k ] = decomp_pred_para( A,b,C_e,d_e,C_i,d_i, i0, v0,lam0,mu0, beta,gamma, tol, rho_uzawa,tol_uzawa )
%% Algorthime de d�compostion par pr�diction en version parall�le
%% pour un probl�me quadratique de type :
%%   minimiser  J(u) := 1/2<A*u,u> - <b,u>
%%      sous     C_e*u == d_e        (e pour contraintes d'�galit�)
%%               C_i*u <= d_i        (i pour contraintes d'in�galit�)
%%
%% Les autres param�tres sont :
%% - i0 : l'unit� � laquelle les contraintes seront affect�es
%% - v0 : une allocation de ressources initiale
%% - lam0 : un prix initial pour les contraintes d'�galit� (mettre 0 si inexistantes)
%% - mu0 : un prix initial pour les contraintes d'in�galit� (mettre 0 si inexistantes) 
%% - beta : appartenant � ]0,1], param�tre pour acc�l�rer la convergence lors de la mise � jour du prix 
%% - gamma : appartenant � ]0,1], param�tre pour acc�l�rer la convergence lors de la mise � jour de l'allocation  
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

V = v0; 
Lam = lam0; % en colonne, les multiplicateurs de Lagrange associ�s aux N sous-probl�mes
Mu = mu0;

%-------------------------------------------------------------------------------------------------------------------------------------------
% Cas 1 : s'il n'y a que des contraintes d'�galit�
if C_i == 0
  % Premi�re it�ration (k = 1)
  
  %  i) Unit� i0
  [U(i0,1),Lam(:,2),~,~] = Uzawa(A(i0,i0),b(i0),C_e(:,i0),V(:,1),0,0, Lam(:,1),0, rho_uzawa,tol_uzawa);
  Lam(:,2) = (1-beta)*Lam(:,1) + beta*Lam(:,2);
  %  ii) Autres unit�s (i != i0)
  for i = [1:(i0-1),(i0+1):N]
      U(i,1) = (b(i) - Lam(:,1)' * C_e(:,i))/A(i,i);
  end
  V(:,2) = (1-gamma)*V(:,1) + gamma*(d_e - (C_e*U(:,1) - C_e(:,i0)*U(i0,1)));
    
  err = tol + 1;

  k = 1;
  while err > tol && k < iter_max
      k = k+1; 
      
      % i) Unit� i0
      [U(i0,k),Lam(:,k+1),~,~] = Uzawa(A(i0,i0),b(i0),C_e(:,i0),V(:,k),0,0, Lam(:,k),0, rho_uzawa,tol_uzawa);
      Lam(:,k+1) = (1-beta)*Lam(:,k) + beta*Lam(:,k+1);
      
      % ii) Autres unit�s (i != i0)
      for i = [1:(i0-1),(i0+1):N]
          U(i,k) = (b(i) - Lam(:,k)' * C_e(:,i))/A(i,i);
      end
      V(:,k+1) = (1-gamma)*V(:,k) + gamma*(d_e - (C_e*U(:,k) - C_e(:,i0)*U(i0,k)));
      
      err = norm(U(:,k) - U(:,k-1));
  end


%-------------------------------------------------------------------------------------------------------------------------------------------
% Cas 2 : s'il n'y a que des contraintes d'in�galit�
elseif C_e == 0 
  % Premi�re it�ration (k = 1)
  
  %  i) Unit� i0
  [U(i0,1),~,Mu(:,2),~] = Uzawa(A(i0,i0),b(i0),0,0,C_i(:,i0),V(:,1), 0,Mu(:,1), rho_uzawa,tol_uzawa);
  Mu(:,2) = (1-beta)*Mu(:,1) + beta*Mu(:,2);
  %  ii) Autres unit�s (i != i0)
  for i = [1:(i0-1),(i0+1):N]
      U(i,1) = (b(i) - Mu(:,1)' * C_i(:,i))/A(i,i);
  end
  V(:,2) = (1-gamma)*V(:,1) + gamma*(d_i - (C_i*U(:,1) - C_i(:,i0)*U(i0,1)));
    
  err = tol + 1;

  k = 1;
  while err > tol && k < iter_max
      k = k+1; 
      
      % i) Unit� i0
      [U(i0,k),~,Mu(:,k+1),~] = Uzawa(A(i0,i0),b(i0),0,0,C_i(:,i0),V(:,k), 0,Mu(:,k), rho_uzawa,tol_uzawa);
      Mu(:,k+1) = (1-beta)*Mu(:,k) + beta*Mu(:,k+1);
      
      % ii) Autres unit�s (i != i0)
      for i = [1:(i0-1),(i0+1):N]
          U(i,k) = (b(i) - Mu(:,k)' * C_i(:,i))/A(i,i);
      end
      V(:,k+1) = (1-gamma)*V(:,k) + gamma*(d_i - (C_i*U(:,k) - C_i(:,i0)*U(i0,k)));
      
      err = norm(U(:,k) - U(:,k-1));
  end

  
%-------------------------------------------------------------------------------------------------------------------------------------------
% Cas 3 : s'il y a des contraintes d'�galit� et d'in�galit�
else
  m = length(d_e) + length(d_i); # nombre de contraintes 
  nb_ega = length(d_e); # nombre de contraintes d'�galite, il y a donc m - nb_ega contraintes d'in�galit�

  % Premi�re it�ration (k = 1)
  
  %  i) Unit� i0
  [U(i0,1),Lam(:,2),Mu(:,2),~] = Uzawa(A(i0,i0),b(i0),C_e(:,i0),V(1:nb_ega,1),C_i(:,i0),V((nb_ega+1):m,1), Lam(:,1),Mu(:,1), rho_uzawa,tol_uzawa);
  Lam(:,2) = (1-beta)*Lam(:,1) + beta*Lam(:,2);
  Mu(:,2) = (1-beta)*Mu(:,1) + beta*Mu(:,2);
  %  ii) Autres unit�s (i != i0)
  for i = [1:(i0-1),(i0+1):N]
      U(i,1) = (b(i) - Lam(:,1)' * C_e(:,i) - Mu(:,1)' * C_i(:,i))/A(i,i);
  end
  V(:,2) = (1-gamma)*V(:,1) + gamma*[d_e - (C_e*U(:,1) - C_e(:,i0)*U(i0,1)) ; d_i - (C_i*U(:,1) - C_i(:,i0)*U(i0,1))];
    
  err = tol + 1;

  k = 1;
  while err > tol && k < iter_max
      k = k+1; 
      
      %  i) Unit� i0
      [U(i0,k),Lam(:,k+1),Mu(:,k+1),~] = Uzawa(A(i0,i0),b(i0),C_e(:,i0),V(1:nb_ega,k),C_i(:,i0),V((nb_ega+1):m,k), Lam(:,k),Mu(:,k), rho_uzawa,tol_uzawa);
      Lam(:,k+1) = (1-beta)*Lam(:,k) + beta*Lam(:,k+1);
      Mu(:,k+1) = (1-beta)*Mu(:,k) + beta*Mu(:,k+1);
      
      %  ii) Autres unit�s (i != i0)
      for i = [1:(i0-1),(i0+1):N]
          U(i,k) = (b(i) - Lam(:,k)' * C_e(:,i) - Mu(:,k)' * C_i(:,i))/A(i,i);
      end
      V(:,k+1) = (1-gamma)*V(:,k) + gamma*[d_e - (C_e*U(:,k) - C_e(:,i0)*U(i0,k)) ; d_i - (C_i*U(:,k) - C_i(:,i0)*U(i0,k))];
      
      err = norm(U(:,k) - U(:,k-1));
  end 

end

if k == iter_max
    disp('ATTENTION : nombre maximal d''it�rations atteint dans ''decomp_pred_para''')
end

endfunction