function [triInds, barycenterCoords] = projectVerticesOnTriangleMesh(X, T, measureLocations)
    
    m = size(measureLocations,1);
    nT = size(T,1);
    nV = size(X,1);
    
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
    
    triInds = zeros(m,1);
    barycenterCoords = zeros(m,3);
    for i =1:m
        point = measureLocations(i,:)';
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
        [~, closestTri] = min(norms(reshape(matrixA*barycens-p,3,[])',2,2));
        barycentersExtracted = reshape(barycens,3,[])';
        triInds(i)=closestTri;
        barycenterCoords(i,:) = barycentersExtracted(closestTri,:);
    end
end