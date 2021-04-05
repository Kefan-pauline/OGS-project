function [u,v,w,it] = alloc(alpha,beta,C,d0,a, b,rho_alloc,itmax,tol_alloc,v0,w0)

%Initialisation
it = 1;
N = size(alpha,1)-1;
nbC = size(C,1);
lambda = zeros(nbC,N);
v = zeros(2,N);;
rho = 0.1;
u = w0;
w = w0;
w_next = w;
d = d0;
tmp = zeros(2,1);
err = 1;
while err > tol_alloc & it < itmax
    
    %Décomposition : resolution des N sous-pbs par l'algo de Uzawa avec inégalité
    for i=1:N
        A = [alpha(i) 0; 0 beta(i)];
        %pour chaque ss pb on resouds avec Uzawa_Inegalite et on stocke les solutions
        [v(:,i),~,lambda(:,i),~] = Uzawa_Inegalite(A,tmp,C,d(:,i),rho,tol_alloc,v0(:,i),itmax,lambda(:,i));
    end
    d = [w(1,1:N)-a;b-w(1,1:N)-w(2,1:N);-b+w(1,1:N)+w(2,1:N)];
    
    %Mise a jour des Allocations
    lambda(1:2,N+1) = [-alpha(N+1)*w(1,N+1) -beta(N+1)*w(2,N+1)]';
    
    for i=1:N+1  
        w_next(:,i) = w(:,i) + rho_alloc*(lambda(1:2,i)-(1/(N+1))*sum(lambda(1:2),2));
    end
    
    
    err = norm(v-v0); 
    it = it+1;
    v0 = v;
   
    % Mise à jour des variables
    w0 = w;
    w = w_next;
    u = w_next;
    u(:,N+1) = -w_next(:,N+1);

end

end

