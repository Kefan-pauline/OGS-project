function [u,v,P,it] = pred(alpha,beta,a,b,Cu,Cv,P,v0, tol_pred,itmax,seq)

N = size(alpha,1)- 1;
u = zeros(2,N);
u0 = u;
v = v0;

%du = zeros(3,N);
diffu = 1; diffv = 1;
it = 0;

while it < itmax && diffu > tol_pred && diffv > tol_pred
    
    % Resolution du sous-probleme N+1
    u(:,N+1) = sum(u,2);
    if seq
        P = [alpha(N+1)*u(1,N+1) beta(N+1)*u(2,N+1)]';
    end
    du = [v(1,:)-a; b-v(1,:)-v(2,:); -b+v(1,:)+v(2,:); a; zeros(1,N); b-a; zeros(1,N)];
    % Resolution des sous-reseaux i = 1,...,N
    for i = 1:N
        A = [alpha(i) 0; 0 beta(i)];
        fun = @(u)Ju(A,P,u,v(:,i));
        u(:,i) = fmincon(fun,u(:,i),Cu,du(:,i));
        
    end
    
    v = v0;
    dv = [u(1,1:N)-a;b-u(1,1:N)-u(2,1:N);-b+u(1,1:N)+u(2,1:N)];
    for i=1:N
        A = [alpha(i) 0; 0 beta(i)];
        fun = @(v)Ju(A,P,u(:,i),v);
        v(:,i)=fmincon(fun,v(:,i),Cv,dv(:,i));
       
    end
    
    if ~seq
        P = [alpha(N+1)*u(1,N+1) beta(N+1)*u(2,N+1)]';
    end   
    diffu = norm(u(:,1:N)-u0,inf);
    u0 = u(:,1:N);
    v0 = v;
    it=it+1;
end
end

