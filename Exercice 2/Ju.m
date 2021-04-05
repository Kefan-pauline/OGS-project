function [J] = Ju(A,p,u,v)
% Fonction objective du probleme de decomposition par les prix
J = 0.5 * v'*A*v+ u'*p;
end

