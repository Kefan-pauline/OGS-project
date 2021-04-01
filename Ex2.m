% Ex 2 : Reseuax de distribution d'eau connectes
clear all;
clc
%% Initialisation
N = 5;
M = 2;
max = 10;
alpha = max*rand(N+1,1);
beta = max*rand(N+1,1);

a = randi(10,1,N);
b = randi(10,1,N)+ a;


C = [-1 0 -1 0 ; 1 1 1 1; -1 -1 -1 -1; 1 0 0 0; -1 0 0 0;0 1 0 0; 0 -1 0 0];
d = [-a;b;-b;a; zeros(1,N); b-a;zeros(1,N)];

%% 1 Decomposition par prix, par quantites et par prediction
% Prix
% Parametre
par.rho = 0.1;
par.ksi = 0.1;
par.tol = 10^-4;
par.itMax = 1000;
par.u0 = ones(4,N+1);
par.p0 = ones(4,1);

[u,p,it] = prix( C, d, alpha, beta, par );
u
%% Quantite
C = [-1 0; 1 1;-1 -1];

w0 = rand(M,N);
w0(:,N+1) = - sum(par.w0,2);

d0 = [w0(1,1:N) - a; b - w0(1,1:N) - w0(2,1:N); -b + w0(1,1:N) + w0(2,1:N)];

par.v0 = ones(M,N);
par.w0 = w0
par.rho = 0.1; 
itmax = 1000; 
tol =10^-3; 


[u,v,w,it] = alloc( C,d0, alpha, beta, a, b, par );

%% Prediction
C = [-1 0 -1 0 ; 1 1 1 1; -1 -1 -1 -1; 1 0 0 0; -1 0 0 0;0 1 0 0; 0 -1 0 0];
par.rho = 0.1;
par.ksi = 0.1;
par.tol = 10^-4;
par.itMax = 1000;
par.u0 = ones(4,N+1);
par.p0 = ones(4,1);






