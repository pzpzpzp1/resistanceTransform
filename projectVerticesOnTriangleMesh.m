function [triInds, barycenterCoords] = projectVerticesOnTriangleMesh(X, T, measureLocations)
    
    m = size(measureLocations,1);
    nT = size(T,1);
    nV = size(X,1);
    for i =1:m
        point = measureLocations(i,:)';
        A33T = reshape(X(reshape(T',[],1),:)',3,3,[]);
        matrixA = sparse(3*nT,3*nT);
        matrixA(1:3*nT*3+3:end) = A33T(1,1,:);
        matrixA(2:3*nT*3+3:end) = A33T(2,1,:);
        matrixA(3:3*nT*3+3:end) = A33T(3,1,:);
        matrixA(3*nT+1:3*nT*3+3:end) = A33T(1,2,:);
        matrixA(3*nT+2:3*nT*3+3:end) = A33T(2,2,:);
        matrixA(3*nT+3:3*nT*3+3:end) = A33T(3,2,:);
        matrixA(2*3*nT+1:3*nT*3+3:end) = A33T(1,3,:);
        matrixA(2*3*nT+2:3*nT*3+3:end) = A33T(2,3,:);
        matrixA(2*3*nT+3:3*nT*3+3:end) = A33T(3,3,:);
        
        p = repmat(point,nT,1);
        
        cvx_begin
            cvx_solver mosek;
            variable barycens(3*nT,1);
            minimize norm(matrixA*barycens-p);
            subject to
                barycens >= 0;
                barycens <= 1;
                sum(reshape(barycens,3,[])) == ones(1,nT);
        cvx_end
        [probzero, closestTri] = min(norms(reshape(matrixA*barycens-p,3,[])',2,2));
        barycentersExtracted = reshape(barycens,3,[])';
        
        norm(X(T(closestTri,:),:)'*barycentersExtracted (closestTri,:)'-point)
        
        
        % project single point onto triangle mesh.
%         point = measureLocations(i,:);
%         triInds = reshape(data.triangles(data.isBoundaryTriangle==1,:)',[],1);
%         bulkMat = permute(reshape(data.vertices(triInds,:)',3,3,[]),[2 1 3]);
%         bulkMat(4,:,:)=1;
%         bulkVec = repmat([point'; 1],1,sum(data.isBoundaryTriangle));
%         barycens = multiSolve(bulkMat, bulkVec);
%         %data.vertices(data.triangles(find(data.isBoundaryTriangle==1,1),:),:)
%         
%         find(norms(permute(multiprod(bulkMat, permute(barycens,[1 3 2])),[1 3 2])-bulkVec,2)<1e-2 & all(barycens>0))
    end
    




end