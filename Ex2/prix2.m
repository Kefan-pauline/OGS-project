function [u,v,P,it] = prix2(alpha,beta,a,b,Cu,Cv,P,v0,rho_prix, tol_prix,itmax)

N = size(alpha,1)- 1;
u = zeros(2,N);
u0 = u;
v = v0;
dv = zeros(3,N);
%du = zeros(3,N);
diffu = 1; diffv = 1;
it = 0;

while it < itmax && diffu > tol_prix && diffv > tol_prix
    du = [v(1,:)-a; b-v(1,:)-v(2,:); -b+v(1,:)+v(2,:); a; zeros(1,N); b-a; zeros(1,N)];
    % Resolution des sous-reseaux i = 1,...,N
    for i = 1:N
        A = [alpha(i) 0; 0 beta(i)];
        fun = @(u)Ju(A,P,u,v(:,i));
        u(:,i) = fmincon(fun,u(:,i),Cu,du(:,i));
        
    end
    % Resolution du sous-probleme N+1
    u(:,N+1) = [alpha(N+1) 0; 0 beta(N+1)]\P
    
    v = v0;
    dv = [u(1,1:N)-a;b-u(1,1:N)-u(2,1:N);-b+u(1,1:N)+u(2,1:N)];
    for i=1:N
        A = [alpha(i) 0; 0 beta(i)];
        fun = @(v)Ju(A,P,u(:,i),v);
        v(:,i)=fmincon(fun,v(:,i),Cv,dv(:,i));
    end
    
    diffu = norm(u(:,1:N)-u0,inf);
    diffv = norm(v-v0,inf);
    u0 = u(:,1:N);
    v0 = v;
    it=it+1;
    
    % Coordination des prix
    P = P + rho_prix * (sum(u(:,1:(end-1)),2) - u(:,end));
end

