function [ KKT1,KKT2,KKT3,KKT4 ] = testKKT( A,b,C,d,u1,lambda1)
% test KKT for the quadratic case with affine inequality constraints
% 1/2<Au,u>-<b,u> such that Cu <= d
% u1 : a matrix of solutions
% lambda1 : a matrix of dual solutions
dim=size(u1,2);
if dim>2
    b=repmat(b,1,dim);
    d=repmat(d,1,dim);
end
% gradient = 0
KKT1=A*u1-b+C'*lambda1;
% respecter la contrainte
KKT2=C*u1-d;
% solution duale positive
KKT3=lambda1;
% ecart complementaire
KKT4=lambda1.*(C*u1-d);

end

