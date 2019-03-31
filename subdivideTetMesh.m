% takes tet mesh in XT, outputs 8 tet subdivision per tet in XoTo.
function [Xo,To] = subdivideTetMesh(X,T)

    nV = size(X,1);
    edges = reshape(T(:,[1 2 1 3 1 4 2 3 3 4 4 2])',2,[])';
    edges = [edges; edges(:,[2 1])];
    edges = edges(edges(:,1)<edges(:,2),:);
    edges = unique(edges,'rows');
    nE = size(edges,1);

    newX = (X(edges(:,1),:)+X(edges(:,2),:))/2;
    vv2newX = sparse(edges(:,1),edges(:,2),1:nE,nV,nV);
    vv2newX=vv2newX+vv2newX';
    vv2newX = vv2newX + nV;

    % eight new tets. 4 on corners. 4 interior split of an octahedron.
    % 1 12 13 14
    % 2 23 21 24
    % 3 31 32 34
    % 4 41 43 42
    % 12 34 13 14
    % 12 34 14 24
    % 12 34 24 23
    % 12 34 23 13

    x1 = T(:,1);
    x2 = T(:,2);
    x3 = T(:,3);
    x4 = T(:,4);
    x12 = vv2newX(sub2ind(size(vv2newX),T(:,1),T(:,2)));
    x13 = vv2newX(sub2ind(size(vv2newX),T(:,1),T(:,3)));
    x14 = vv2newX(sub2ind(size(vv2newX),T(:,1),T(:,4)));
    x23 = vv2newX(sub2ind(size(vv2newX),T(:,2),T(:,3)));
    x34 = vv2newX(sub2ind(size(vv2newX),T(:,3),T(:,4)));
    x24 = vv2newX(sub2ind(size(vv2newX),T(:,4),T(:,2)));

    To = [x1 x12 x13 x14...
    x2 x23 x12 x24...
    x3 x13 x23 x34...
    x4 x14 x34 x24...
    x12 x34 x13 x14...
    x12 x34 x14 x24...
    x12 x34 x24 x23...
    x12 x34 x23 x13...
    ]';
    To = reshape(To,4,[])';

    Xo = [X; newX];

end

% displays 8 consecutive tets which form one original tet.
%         figure; hold all; axis off; axis equal; rotate3d on;
%         start = randi(size(T,1));
%         for i = 1:8
%             tet = To(start*8+i,:);
%             tris = reshape(tet([1 2 3 1 3 4 1 4 2 2 4 3]),3,[])';
%             patch('Faces',tris,'Vertices',Xo,'FaceColor','green','EdgeColor','black');
%             pause
%         end
        