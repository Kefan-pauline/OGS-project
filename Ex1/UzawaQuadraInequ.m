function [ u0,miu0,it ] = UzawaQuadraInequ( A,b,C,d,u0,miu0,rho,itMax,tol )
% Uzawa algorithm for the quadratic case with affine inequality constraints
% 1/2<Au,u>-<b,u> such that Cu <= d
%   u0 : initial point
%   miu0 : initial price
%   rho : pas de gradient
it=1; 
err=1;
while it<itMax && err>tol 
    u1=A\(b-C'*miu0); 
    miu0= max(0,miu0+rho*(C*u1-d));
    err=norm(u1-u0);
    u0=u1;
    it=it+1; 
end
%fprintf('Uzawa en %d iterations \n',it);
end

