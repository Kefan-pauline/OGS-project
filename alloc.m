function [ u, v, w, it ] = alloc( C,d0,alpha,beta,a, b, par )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% Parametres
v0 = par.v0;
w0 = par.w0;
rho = par.rho;
itMax = par.itMax;
tol = par.tol;

%Initialisation
it=1;
err=1;


N = size(alpha,1)-1;
M = size(C,1);


v = v0;
u = w0;
w = w0;
w_next = w;

d = d0;
tmp = zeros(2,1);
mu0 = zeros(M,N+1);
while (err > tol && it < itMax)
    
    %Décomposition : resolution des N sous-pbs
    for i=1:N
        A = [alpha(i) 0; 0 beta(i)];
        [v(:,i),mu(:,i),~] = Uzawa(A,tmp,0,0,C,d(:,i),0,mu0(:,i),rho,tol);
    end
    
    d = [w(1,1:N)-a;b-w(1,1:N)-w(2,1:N);-b+w(1,1:N)+w(2,1:N)];
    
    mu(1:2,N+1) = [-alpha(N+1)*w(1,N+1);  -beta(N+1)*w(2,N+1)];
    
    for i=1:N
        w_next(:,i) = w(:,i) + rho.*(mu(1:2,i)-(1/N).*sum(mu(1:2,:),2));
    end
    
    it = it +1;
    err = norm(w_next-w);
    
    % mise a jour des variables
    w = w_next;
    u = w_next;
    u(:,N+1) = -w_next(:,N+1);
end

end

