%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Exercice 2 Kefan S., Hiba S., Victor K. Vinh N.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialisation
%close all;
%clear all;
maxi = 10;
N = 10;
alpha = 0 + (maxi-0)*rand(N+1,1);
beta = 0 + (maxi-0)*rand(N+1,1);
a = randi(maxi,1,N);
b = randi(maxi,1,N) + a;

%% Décomposition par le prix 
Cu = [-1 0;1 1;-1 -1 ; 1 0;-1 0 ;0 1; 0 -1];
v0 = ones(2,N);


Cv = [-1 0; 1 1; -1 -1];
rho_prix =0.01;        % pas de l'algo de descente de gradient
nitermax_prix =1000; % nb max iterations prix
tol_prix =10^-3;       % tolerance prix

P0 = ones(2,1);

[u,v,P,it] = prix2(alpha,beta,a,b,Cu,Cv,P0,v0,rho_prix,tol_prix,nitermax_prix);
it
u

%% Question 3:
% C'est long a executer pour cette section
maxi = 10;
inter = 20:10:200;
temps = inter;
iter = inter;

Cu = [-1 0;1 1;-1 -1 ; 1 0;-1 0 ;0 1; 0 -1];
Cv = [-1 0; 1 1; -1 -1];
rho_prix = [1,0.1,0.05,0.01,0.001];         
nitermax_prix = 1000;                       %nb max iterations prix
tol_prix = [10^-3,10^-3,10^-3,10^-3,10^-3]; %precision prix
P0 = ones(2,1);
j=1;
for N = inter
    alpha = 0 + (maxi-0)*rand(N+1,1);
    beta = 0 + (maxi-0)*rand(N+1,1);
    a = randi(maxi,1,N);
    r = randi(maxi,1,N);
    b = r + a;
    v0 = ones(2,N);
    for i=1:N
        v0(1,i) = randi([a(i) maxi],1);
        v0(2,i) = b(i) - v0(1,i);
    end
    
    %1
    tic;
    [u,v,P,it] = prix2(alpha,beta,a,b,Cu,Cv,P0,v0,rho_prix(1),tol_prix(1),nitermax_prix);
    temps(j) = toc;
    iter(j) = it;
    %2
    tic;
    [u,v,P,it] = prix2(alpha,beta,a,b,Cu,Cv,P0,v0,rho_prix(2),tol_prix(2),nitermax_prix);
    temps2(j) = toc;
    iter2(j) = it;
    %3
    tic;
    [u,v,P,it] = prix2(alpha,beta,a,b,Cu,Cv,P0,v0,rho_prix(3),tol_prix(3),nitermax_prix);
    temps3(j) = toc;
    iter3(j) = it;
    %4
    tic;
    [u,v,P,it] = prix2(alpha,beta,a,b,Cu,Cv,P0,v0,rho_prix(4),tol_prix(4),nitermax_prix);
    temps4(j) = toc;
    iter4(j) = it;
    %5
    tic;
    [u,v,P,it] = prix2(alpha,beta,a,b,Cu,Cv,P0,v0,rho_prix(5),tol_prix(5),nitermax_prix);
    temps5(j) = toc;
    iter5(j) = it;
    j=j+1;
end
%%
figure(1)
plot(inter, temps, inter, temps2, inter, temps3,inter, temps4,inter, temps5)
xlabel('N')
ylabel('Temps (s)')
title('Temps en fonction de N')
legend('\rho = 1, \epsilon = 10^{-3}', '\rho = 0.1, \epsilon = 10^{-3}','\rho = 0.05, \epsilon = 10^{-3}','\rho = 0.01, \epsilon = 10^{-3}','\rho = 0.001, \epsilon = 10^{-3}')
%%
% figure(2)
% plot(inter, iter,inter, iter2,inter, iter3,inter, iter4,inter, iter5)
% title('Nombre d iteration en fonction de N')
% xlabel('N')
% ylabel('Nombre d iteration')
% legend('\rho = 1, \epsilon = 10^{-3}', '\rho = 0.1, \epsilon = 10^{-3}','\rho = 0.05, \epsilon = 10^{-3}','\rho = 0.01, \epsilon = 10^{-3}','\rho = 0.001, \epsilon = 10^{-3}')
%% allocation
%close all;
%clear;
N = 5;
alpha = 0 + (10-0)*rand(N+1,1);
beta = 0 + (10-0)*rand(N+1,1);
C = [-1 0; 1 1;-1 -1];
a = randi(10,1,N);
b = randi(10,1,N) + a;

w0 = randi(10,2,N);
w0(:,N+1) = -sum(w0,2);
d0 = [w0(1,1:N)-a;b-w0(1,1:N)-w0(2,1:N);-b+w0(1,1:N)+w0(2,1:N)];

rho_alloc = 0.1;      %pas de l'algo de descente de gradient
itmax_alloc = 1000;   %nb max iterations allocation
tol_alloc = 10^-3;    %precision allocation
v0 = ones(2,N);

[u,v,w,it] = alloc(alpha,beta,C,d0,a, b,rho_alloc,itmax_alloc,tol_alloc,v0,w0);
u
%% Prediction
N = 5;
alpha = 0 + (10-0)*rand(N+1,1);
beta = 0 + (10-0)*rand(N+1,1);
a = randi(10,1,N);
b = randi(10,1,N) + a;
Cu = [-1 0;1 1;-1 -1 ; 1 0;-1 0 ;0 1; 0 -1];
v0 = ones(2,N);
% Variable pour utiliser la version sequentielle ou non
seq = true;

Cv = [-1 0; 1 1; -1 -1];
nitermax_pred = 1000;   % nb max iterations prix
tol_prix =10^-3;        % tolerance prix

P0 = ones(2,1);

[u,v,P,it] = pred(alpha,beta,a,b,Cu,Cv,P0,v0,tol_prix,nitermax_pred,seq);
it
u
