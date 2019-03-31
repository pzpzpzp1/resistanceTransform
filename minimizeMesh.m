% Given a removes unreferenced vertices of X and re-indexes T accordingly. 
% ut is the index to transfer a scalar field from X to a scalar field on
% Xo: phi_o = phi(ut);
function [Xo,To,ut]=minimizeMesh(X,T)
    
    nV = size(X,1);
    nT = size(T,1);
    
    [ut, ia, ic] = unique(T);
    Xo = X(ut,:);
    nXo=size(Xo,1);
    indexXo = 1:nXo;
    To = reshape(indexXo(ic),[],3);
    
end