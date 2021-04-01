function [ u, p, it ] = prix( C, d,alpha, beta, par )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
it = 1;
err = 1;
n = size(alpha,1); % nombre de sous-reseaux
m = size(C,1); % number of constraints

u0 = par.u0;
p0 = par.p0;

u = u0;
p = p0;

while (err > par.tol && it < par.itMax)
    for i = 1:n
        A = diag([0 0 alpha(i) beta(i)]);
        if i ~= n % Résolution des sous-problèmes i = 1,...,N
            [u(:,i),~,it] = ArrowHurwiczQuadra( A,-p,C,d(:,i),u(:,i),par.ksi,par.rho,par.itMax,par.tol );
        else % Résolution de sous-problèmes i = N+1
            u(3:4,i) = 0;
            u(1:2,i) = A(3:4,3:4)\p(3:4);  
        end
    end
    it = it + 1;
    
    % Coordination
    p = p + par.rho * (sum(u(:,1:end-1),2)-u(:,end));
    
    err = norm(u - u0); 
    u0 = u;
end

end

