% data is tet mesh data
% perform FEM sim of injecting current into sourcePos, and out of sinkPos
% and read off resulting voltages.
function measuredVoltages = measureVoltages(data, tetra, Laplacian, sourcePos, sinkPos, current, measureLocations, resistivity)
    %{
    % extract boundary tet mesh.
    bT = data.triangles(data.isBoundaryTriangle==1,:);
    bX = data.vertices;
    [bX,bT]=minimizeMesh(bX,bT);
    FV.faces = bT;
    FV.vertices = bX;
    [distances,surface_points] = point2trimesh(FV, 'QueryPoints', measureLocations);
    if debugging
        figure; hold all; rotate3d on; 
        ptc=patch('Faces',FV.faces,'Vertices',FV.vertices,'FaceColor','green','EdgeColor','black');
        alpha(ptc,.1);
%         scatter3(measureLocations(:,1),measureLocations(:,2),measureLocations(:,3),'r.'); 
        scatter3(measureLocations(23,1),measureLocations(23,2),measureLocations(23,3),'r.');
        scatter3(surface_points(23,1),surface_points(23,2),surface_points(23,3),'b.');
%         scatter3(data.vertices(data.isBoundaryVertex==1,1),data.vertices(data.isBoundaryVertex==1,2),data.vertices(data.isBoundaryVertex==1,3),'g');
    end
    %}

    % locate source and sink on tet mesh.
    [sourceTet, barycenterSourceCoords ] = pointLocation(tetra,sourcePos);
    [sinkTet, barycenterSinkCoords] = pointLocation(tetra,sinkPos);

    % check that source/sink points are on the boundary of some tet. 
    assert(numel(find(abs(barycenterSourceCoords)<1e-6))==1);
    assert(numel(find(abs(barycenterSinkCoords)<1e-6))==1);
    
    sourceSupportIndicator = (abs(barycenterSourceCoords) > 1e-6);
    sinkSupportIndicator = (abs(barycenterSinkCoords) > 1e-6);
    sourceSupport = data.tetrahedra(sourceTet,sourceSupportIndicator);
    sinkSupport = data.tetrahedra(sinkTet,sinkSupportIndicator);

    % construct current vector by distributing source and sink
    sourceDists = norms(sourcePos - data.vertices(sourceSupport,:),2,2);
    sinkDists = norms(sinkPos - data.vertices(sinkSupport,:),2,2);
    
    sWeights = (1./sourceDists)/sum(1./sourceDists);
    tWeights = (1./sinkDists)/sum(1./sinkDists);
    
    injectedCurrent = zeros(data.numVertices,1);
    injectedCurrent(sourceSupport) = injectedCurrent(sourceSupport) + current * sWeights;
    injectedCurrent(sinkSupport) = injectedCurrent(sinkSupport) - current * tWeights;
    
    % compute voltages on tet vertices.
    inds = 2:numel(injectedCurrent);
    phiInds = [(Laplacian(inds,inds) \ injectedCurrent(inds))*resistivity]; % base voltage 0 at first vertex.
    phi = zeros(numel(injectedCurrent),1);
    phi(inds)=phiInds;
    phi = phi - min(phi);
    
    if false
        % debug laplacian
%         ind = 2;
%         [V,D] = eigs(Laplacian,ind,'smallestabs'); V = V - min(V); V = V./max(V);
%         scatter3(data.vertices(:,1),data.vertices(:,2),data.vertices(:,3),200,[V(:,ind) V(:,ind) V(:,ind)]);
        
        % plot results
        colors = phi; colors = colors - min(colors); colors = colors/max(colors);
        currentColors = injectedCurrent==0;
        figure; hold all; rotate3d on;
%         ptc=patch('Faces',data.triangles(data.isBoundaryTriangle==1,:),'Vertices',data.vertices,'FaceColor','green','EdgeColor','none'); alpha(ptc,.1);
        scatter3(sourcePos(:,1),sourcePos(:,2),sourcePos(:,3),200,'g','filled');
        scatter3(sinkPos(:,1),sinkPos(:,2),sinkPos(:,3),200,'r','filled');
        scatter3(data.vertices(:,1),data.vertices(:,2),data.vertices(:,3),200,[colors colors colors]);
%         scatter3(data.vertices(:,1),data.vertices(:,2),data.vertices(:,3),200,[currentColors currentColors currentColors]);
%         patch('Faces',data.triangles(data.tetsToTriangles(sourceTet,:),:),'Vertices',data.vertices,'FaceColor','blue','EdgeColor','black'); alpha(ptc,.1);
%         patch('Faces',data.triangles(data.tetsToTriangles(sinkTet,:),:),'Vertices',data.vertices,'FaceColor','blue','EdgeColor','black'); alpha(ptc,.1);
    end
    
    % interpolate voltages linearly to sample at measured locations
    measuredVoltages = SampleScalarFieldFromTetMesh(data, tetra, phi, measureLocations, 0);
end