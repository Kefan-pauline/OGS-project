%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Exercice 1 (Kefan S., Hiba S., Victor K. Vinh N.  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

% exo 1
% 1.
[U1,U2] = meshgrid(-5:0.01:5,-5:0.01:5);
LM = U1 + 2*U2 <= 0 & U2 <=0;
U1(~LM) = NaN;
U2(~LM) = NaN;
J = (U1.^2+U2.^2)/2-U1+U2;

figure
surf(U1, U2, J,'FaceAlpha',0.5,'EdgeColor','none');
view(2);
xlabel('u1');
ylabel('u2');
title('Les solutions admissibles');

figure
surf(U1, U2, J,'FaceAlpha',0.5,'EdgeColor','none');
xlabel('u1');
ylabel('u2');
zlabel('J');
title('La fonction');
%%
% 2. 
% Set up for N = 4
N = 4;
A = sparse(1:N,1:N,ones(1,N),N,N);
b = repmat([1;-1],N/2,1);
C = A + sparse(2:N,1:N-1,2*ones(1,N-1),N,N)';
d = zeros(N,1);

% Uzawa
u0 = zeros(N,1);
miu0 = zeros(N,1);
rho_U = .2;
itMax = 1000;
tol = 10^(-4);
[ u,miu,it ] = UzawaQuadraInequ( A,b,C,d,u0,miu0,rho_U,itMax,tol );
%%
% 2.
% Decomposition par prix
u0 = zeros(N,1);
miu0 = zeros(N,1);
rho_U = .1;
itMax = 1000;
tol = 10^(-4);
[ u_p,p_p,it_p,KKT_p ] = DecompPrice( A,b,C,d,u0,miu0,rho_U,itMax,tol );
%%
% 2.
% Decomposition par quantite
u0 = zeros(N,1);
v0 = zeros(N,N);
rho = .1;
rho_U = .1;
itMax = 1000;
tol = 10^(-4);
tol_U = 10^(-5);
[ u_q,p_q,it_q,KKT_q ] = DecompQuant( A,b,C,d,u0,v0,rho,rho_U,tol_U,itMax,tol );
%%
% 2.
% Decomposition par prediction (parallele)
u0 = zeros(N,1);
i0 = N-1;
v0 = zeros(N,1);
p0 = zeros(N,1);
gamma = .5;
beta = .5;
rho_U = 0.1;
itMax = 1000;
tol = 10^(-4);
tol_U = 10^(-5);
[ u_pred,p_pred,it_pred,KKT_pred ] = DecompPredPara(A,b,C,d,i0,u0,v0,p0,gamma,beta,rho_U,tol_U,itMax,tol );
%%
% 3. 
% KKT decomposition par prix
[ KKT1,KKT2,KKT3,KKT4 ]=testKKT(A,b,C,d,u_p,p_p);
%%
% 3.
% KKT decomposition par quantite
[ KKT1,KKT2,KKT3,KKT4 ]=testKKT(A,b,C,d,u_q,reshape(p_q(:,1,:),[N,it_q]));
%%
% 3.
% KKT decomposition par predction
[ KKT1,KKT2,KKT3,KKT4 ]=testKKT(A,b,C,d,u_pred,p_pred);
%%
% 3.
% Comparaison des resultats
N=4;
tol=10^-4;
tol_U=0.1*tol;
rho=.5;
rho_U=.1;
[ It,Temps] = Comparaison( N,tol,tol_U,rho,rho_U,1 );
%%
% 4.
% Etude de la complexite (N=2,4,...,26)
% Ce bout de code met environ 20s a executer
tol=10^-4;
tol_U=0.1*tol;
rho=.2;
rho_U=.1;
list_N=zeros(1,13);
list_it=zeros(3,13);
list_temps=zeros(3,13);
for k=1:13
    list_N(k)=2*k;
    [ list_it(:,k),list_temps(:,k)] = Comparaison( list_N(k),tol,tol_U,rho,rho_U,0 );
end
figure
loglog(list_N,list_it)
title('It en fonction du N')
xlabel('N')
ylabel('It')
legend('Par prix','Par quantité','Par prediction')
figure
loglog(list_N,list_temps)
title('Temps en fonction du N ')
xlabel('N')
ylabel('Temps')
legend('Par prix','Par quantité','Par prediction')
%%
% 5.
% Resolution pour N = 200
% Ce bout de code met environ 13s a executer
N=200;
rho=0.1; 
rho_U=0.1;
tol=10^(-4);
tol_U=10^(-5);
[ it,temps ] = Comparaison( N,tol,tol_U,rho,rho_U,0 );
%%
% 5.
% Resolution pour N = 251
% Ce bout de code met environ 30s a executer
N=251;
rho=0.1; 
rho_U=0.1;
tol=10^(-4);
tol_U=10^(-5);
[ it,temps ] = Comparaison( N,tol,tol_U,rho,rho_U,0 );
%%
% 5.
% N = 200 en fonction du tol
% Ce bout de code met environ 40s a executer
N=200;
rho=0.1; 
rho_U=0.1;
list_tol=zeros(1,4);
list_it=zeros(3,4);
list_temps=zeros(3,4);
tol_U=10^(-5);
for k=2:5
    list_tol(k-1)=10^(-k);
    [ list_it(:,k-1),list_temps(:,k-1) ] = Comparaison( N,list_tol(k-1),tol_U,rho,rho_U,0 );
end
figure
loglog(list_tol,list_it)
title(['It en fonction du tol N= ' num2str(N)])
xlabel('tol')
ylabel('It')
legend('Par prix','Par quantité','Par prediction')
figure
loglog(list_tol,list_temps)
title(['Temps en fonction du tol N= ' num2str(N)])
xlabel('tol')
ylabel('Temps')
legend('Par prix','Par quantité','Par prediction')
%%
% 5.
% N = 251 en fonction du tol
% Ce bout de code met environ 57s a executer
N=251;
rho=0.08; 
rho_U=0.1;
list_tol=zeros(1,4);
list_it=zeros(3,4);
list_temps=zeros(3,4);
tol_U=10^(-5);
for k=2:5
    list_tol(k-1)=10^(-k);
    [ list_it(:,k-1),list_temps(:,k-1) ] = Comparaison( N,list_tol(k-1),tol_U,rho,rho_U,0 );
end
figure
loglog(list_tol,list_it)
title(['It en fonction du tol N= ' num2str(N)])
xlabel('tol')
ylabel('It')
legend('Par prix','Par quantité','Par prediction')
figure
loglog(list_tol,list_temps)
title(['Temps en fonction du tol N= ' num2str(N)])
xlabel('tol')
ylabel('Temps')
legend('Par prix','Par quantité','Par prediction')
%%
% 5.
% N = 8 en fonction du rho_U
% Ce bout de code met environ 10s a executer
N=8;
rho = 0.1;
tol = 10^(-4);
tol_U = 0.1*tol;
list_rhoU=[0.01 0.03 0.05 0.07 0.09 0.1];
list_it=zeros(3,6);
list_temps=zeros(3,6);
for k=1:6
    [ list_it(:,k),list_temps(:,k) ] = Comparaison( N,tol,tol_U,rho,list_rhoU(k),0 );
end
figure
loglog(list_rhoU,list_it)
title(['It en fonction du rho_U N= ' num2str(N)])
xlabel('rho_U')
ylabel('It')
legend('Par prix','Par quantité','Par prediction')
figure
loglog(list_rhoU,list_temps)
title(['Temps en fonction du rho_U N= ' num2str(N)])
xlabel('rho_U')
ylabel('Temps')
legend('Par prix','Par quantité','Par prediction')
%%
% 5.
% N = 9 en fonction du rho_U
% Ce bout de code met environ 13s a executer
N=9;
rho = 0.1;
tol = 10^(-4);
tol_U = 0.1*tol;
list_rhoU=[0.01 0.03 0.05 0.07 0.09 0.1];
list_it=zeros(3,6);
list_temps=zeros(3,6);
for k=1:6
    [ list_it(:,k),list_temps(:,k) ] = Comparaison( N,tol,tol_U,rho,list_rhoU(k),0 );
end
figure
loglog(list_rhoU,list_it)
title(['It en fonction du rho_U N= ' num2str(N)])
xlabel('rho_U')
ylabel('It')
legend('Par prix','Par quantité','Par prediction')
figure
loglog(list_rhoU,list_temps)
title(['Temps en fonction du rho_U N= ' num2str(N)])
xlabel('rho_U')
ylabel('Temps')
legend('Par prix','Par quantité','Par prediction')
%%
% 5.
% N = 8 en fonction du rho
N=8;
rho_U = 0.1;
tol = 10^(-4);
tol_U = 0.1*tol;
list_rho=[0.01 0.03 0.05 0.07 0.09 0.1];
list_it=zeros(3,6);
list_temps=zeros(3,6);
for k=1:6
    [ list_it(:,k),list_temps(:,k) ] = Comparaison( N,tol,tol_U,list_rho(k),rho_U,0 );
end
figure
plot(list_rhoU,list_it)
title(['It en fonction du rho N= ' num2str(N)])
xlabel('rho')
ylabel('It')
legend('Par prix','Par quantité','Par prediction')
figure
plot(list_rhoU,list_temps)
title(['Temps en fonction du rho N= ' num2str(N)])
xlabel('rho')
ylabel('Temps')
legend('Par prix','Par quantité','Par prediction')
%%
% 5.
% N = 9 en fonction du rho
N=9;
rho_U = 0.1;
tol = 10^(-4);
tol_U = 0.1*tol;
list_rho=[0.01 0.03 0.05 0.07 0.09 0.1];
list_it=zeros(3,6);
list_temps=zeros(3,6);
for k=1:6
    [ list_it(:,k),list_temps(:,k) ] = Comparaison( N,tol,tol_U,list_rho(k),rho_U,0 );
end
figure
plot(list_rhoU,list_it)
title(['It en fonction du rho N= ' num2str(N)])
xlabel('rho')
ylabel('It')
legend('Par prix','Par quantité','Par prediction')
figure
plot(list_rhoU,list_temps)
title(['Temps en fonction du rho N= ' num2str(N)])
xlabel('rho')
ylabel('Temps')
legend('Par prix','Par quantité','Par prediction')
%%
% 6.
% Ajout d'un sur-diagonal
N = 4;
alpha = -2;
A = sparse(1:N,1:N,ones(1,N),N,N)+ sparse(2:N,1:N-1,-alpha*ones(1,N-1),N,N)';
b = repmat([1;-1],N/2,1);
C = A + sparse(2:N,1:N-1,2*ones(1,N-1),N,N)';
d = zeros(N,1);
% Uzawa
u0 = ones(N,1);
miu0 = ones(N,1);
rho_U = .2;
itMax = 1000;
tol = 10^(-7);
[ u,miu,it ] = UzawaQuadraInequ( A,b,C,d,u0,miu0,rho_U,itMax,tol );
% Test KKT
[ KKT1,KKT2,KKT3,KKT4 ]=testKKT(A,b,C,d,u,miu);