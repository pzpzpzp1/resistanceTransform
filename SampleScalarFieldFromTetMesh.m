%% given tetmesh (data,tetra) and scalar field phi on vertices, measure samples of phi at measureLocations.
%% if measureLocations isn't inside the mesh to within floating point accuracy, return Nan as measured data.
% all means compute interpolated samples for all measure locations, even
% those that fail pointlocation. costs more time.
function measuredVoltages = SampleScalarFieldFromTetMesh(data, tetra, phi, measureLocations, all);

[measureTetInds, BarycenterCoords] = pointLocation(tetra,measureLocations);
measureTetIndsCleaned = measureTetInds;
failLocations = isnan(measureTetInds);
measureTetIndsCleaned(isnan(measureTetInds))=1;
data.tetrahedra(measureTetIndsCleaned,:)

measuredVoltages = sum(BarycenterCoords.*phi(data.tetrahedra(measureTetIndsCleaned,:)),2);

% figure; hold all; rotate3d on; 
% scatter3(measureLocations(2157,1),measureLocations(2157,2),measureLocations(2157,3),'r.')
% ptc=patch('Faces',data.triangles(data.isBoundaryTriangle==1,:),'Vertices',data.vertices,'FaceColor','green','EdgeColor','black');
% alpha(ptc,.1);

% barycenterSourceCoords = data.vertices(data.tetrahedra(measureTetInds,:),:)'\sourcePos';

if all
    [Xo,To,reindex]=minimizeMesh(data.vertices,data.triangles(data.isBoundaryTriangle==1,:));
    restrictedPhi = phi(reindex);
    [triangleInds, failedBarycenterCoords] = projectVerticesOnTriangleMesh(Xo,To,measureLocations(failLocations,:));
    measuredVoltages(failLocations) = sum(failedBarycenterCoords.*restrictedPhi(To(triangleInds,:)),2);
end

end