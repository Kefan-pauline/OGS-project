
function [X,obj,lambda,niter] = Uzawa_Inegalite(A,b,C,d,rho,epsilon,X0,nitermax,lambda0)
% Méthode de Uzawa :
% - A, b, C,d : données du problème
% - u0 : point initial
% - lambda0 : le multiplicateur de lagrange  initial
% - eps : précision voulue
% - nitermax : le nombre d'itérations max
% - rho : le pas
niter=0;
X=X0;
lambda=lambda0;
X = X0;
X_next = X0+epsilon+1; % pour entrer dans le while

while (niter < nitermax && norm(X-X_next,1)>epsilon)
    X_next = X;
    
    X = pinv(A)*(b-C'*lambda); %on resouds le pb global sur les u
    lambda0 = max(0,lambda + rho*(C*X-d)); %descente de gradient sur les lambda
    lambda = lambda0;
    niter = niter+1;
end
%b
obj = 0.5*X'*A*X-X'*b ;%J(A,b,u);
%niter
end
