% Given a removes unreferenced vertices of X and re-indexes T accordingly. 
% idx is the index to transfer a scalar field from X to a scalar field on
% Xo: phi_o = phi(idx);
function [Xo,To,idx]=minimizeHexMesh(X,T)
    
    nV = size(X,1);
    nT = size(T,1);
    
    [idx, ia, ic] = unique(T);
    Xo = X(idx,:);
    nXo=size(Xo,1);
    indexXo = 1:nXo;
    To = reshape(indexXo(ic),[],8);
    
end