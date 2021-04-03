function [ It,Temps ] = Comparaison( N,tol,tol_U,rho,rho_U,affichage )
% Comparer les 3 algorithmes 
%   N : dimension
%   tol : tolerance pour les algorithmes de decomposition
%   tol_U : tolerance pour Uzawa
%   rho : le pas pour la decomposition par quantite
%   rho_U : le pas pour Uzawa
%   affichage = 1 pour voir le graphe qui affiche les solutions de chaque
%       iteration, les points sont en rouge si les condition KKT non verifiees,
%       en vert sinon

% Initialisation
A = sparse(1:N,1:N,ones(1,N),N,N);
C = A + sparse(2:N,1:N-1,2*ones(1,N-1),N,N)';
d = zeros(N,1);
if mod(N,2)
    b = [repmat([1;-1],(N-1)/2,1); 1];
else
    b = repmat([1;-1],N/2,1);
end
Temps = zeros(3,1);

% Initialisation des parametres
u0=zeros(N,1);
p0=zeros(N,1);
v0_q=zeros(N,N);
itMax=1000;
i0=1;
v0=zeros(N,1);

tic
[ u_p,p_p,it_p,KKT_p ] = DecompPrice( A,b,C,d,u0,p0,rho_U,itMax,tol );
Temps(1)=toc;

tic
[ u_q,p_q,it_q,KKT_q ] = DecompQuant( A,b,C,d,u0,v0_q,rho,rho_U,tol_U,itMax,tol );
Temps(2)=toc;

tic
[ u_pred,p_pred,it_pred,KKT_pred] = DecompPredPara(A,b,C,d,i0,u0,v0,p0,.5,.5,rho_U,tol_U,itMax,tol );
Temps(3)=toc;

It = [it_p it_q it_pred];

if(affichage==1)
    % exact solution
    [ u,~,~ ] = UzawaQuadraInequ( A,b,C,d,zeros(N,1),zeros(N,1),.2,itMax,tol );
    % prix
    figure
    myColors = zeros(it_p-1, 3); % List of rgb colors for every data point
    rowsToSetGreen = KKT_p(2:it_p) == 1;
    rowsToSetRed = KKT_p(2:it_p)==0;
    myColors(rowsToSetGreen, :) = repmat([0,1,0],sum(rowsToSetGreen),1);
    myColors(rowsToSetRed, :) = repmat([1,0,0],sum(rowsToSetRed),1);
    scatter(2:it_p,sum((u_p(:,2:it_p)-repmat(u,1,it_p-1)).^2),10,myColors,'filled');
    title(['Décomposition par prix N=' num2str(N)]);
    xlabel('it');
    ylabel('||u-u*||^2');
    legend('KKT non vérifié');
    % quantite
    figure
    myColors = zeros(it_q-1, 3); % List of rgb colors for every data point
    rowsToSetGreen = KKT_q(2:it_q) == 1;
    rowsToSetRed = KKT_q(2:it_q)==0;
    myColors(rowsToSetGreen, :) = repmat([0,1,0],sum(rowsToSetGreen),1);
    myColors(rowsToSetRed, :) = repmat([1,0,0],sum(rowsToSetRed),1);
    scatter(2:it_q,sum((u_q(:,2:it_q)-repmat(u,1,it_q-1)).^2),10,myColors,'filled');
    title(['Décomposition par quantité N=' num2str(N)]);
    xlabel('it');
    ylabel('||u-u*||^2');
    legend('KKT non vérifié');   
    % prediction
    figure
    myColors = zeros(it_pred-1, 3); % List of rgb colors for every data point
    rowsToSetGreen = KKT_pred(2:it_pred) == 1;
    rowsToSetRed = KKT_pred(2:it_pred)==0;
    myColors(rowsToSetGreen, :) = repmat([0,1,0],sum(rowsToSetGreen),1);
    myColors(rowsToSetRed, :) = repmat([1,0,0],sum(rowsToSetRed),1);
    scatter(2:it_pred,sum((u_pred(:,2:it_pred)-repmat(u,1,it_pred-1)).^2),10,myColors,'filled');
    title(['Décomposition par prediction N=' num2str(N)]);
    xlabel('it');
    ylabel('||u-u*||^2');
    legend('KKT non vérifié');
end

end

