
%% declare misc parameters
resistivity = 1;
I = 1;
debugging = 1;
premade = 1;

%% load tet mesh
if premade
    files = dir('meshes/vtkComplementTetMeshes');
    filename = files(randi(numel(files)-2)+2).name;
    [ext,fname,ext] = fileparts(filename);
    fnamec = [fname '_c'];
else
    files = dir('../../10k_tetmesh/');
    filename = files(randi(numel(files)-2)+2).name;
    [ext,fname,ext] = fileparts(filename);
    fnamec = [fname '_c'];

    generateComplementMesh(fname);
end

%% LOAD COMPLEMENT TET MESH
fid = fopen(['meshes/vtkTetMeshes/' fnamec '.vtk']);
[CVerts, CTets] = loadVTKTET(fid);
fclose(fid);

tetra = triangulation(CTets, CVerts);
[TetL, TetM] = GeometricPrimalLM(tetra); % laplacian and mass matrix
Ctetdata = getTetDataRT(CTets,CVerts,1);

%% Solve poisson problem for electrical properties
% isolate rectangle boundary
rectK = convhull(CVerts);
rectVert = rectK(1);
boundaryGraph = graph(Ctetdata.BoundaryEdges(:,1),Ctetdata.BoundaryEdges(:,2));
bins = conncomp(boundaryGraph); 
rectBin = bins(rectVert);
rectVerts = find(bins==1);

% random source/sink selection
twoRandElems = randsample(rectVerts,2,false);
randomSource = twoRandElems(1);
randomSink = twoRandElems(2);

if debugging
    figure; hold all; axis equal; rotate3d on;
    scatter3(Ctetdata.vertices(rectVerts,1),Ctetdata.vertices(rectVerts,2),Ctetdata.vertices(rectVerts,3),100,'k.');
    scatter3(Ctetdata.vertices(randomSource,1),Ctetdata.vertices(randomSource,2),Ctetdata.vertices(randomSource,3),500,'g.');
    scatter3(Ctetdata.vertices(randomSink,1),Ctetdata.vertices(randomSink,2),Ctetdata.vertices(randomSink,3),500,'r.');
end

% div sigma grad phi = I delta_source - I delta_sink. voltages are phi.
currentInjection = zeros(Ctetdata.numVertices,1);
currentInjection(randomSource)=I;
currentInjection(randomSink)=-I;
inds = 1:numel(currentInjection); inds(randomSink)=[];
phi = [(TetL(inds,inds) \ currentInjection(inds))*resistivity]; % base voltage 0 at first vertex.
phi(randomSink:(end+1)) = [0; phi(randomSink:end)];

if debugging
    colors = phi;
    colors = colors-min(colors);
    colors = colors/max(colors);
    colors = [colors colors colors];
    figure; hold all; axis equal; rotate3d on;
    scatter3(Ctetdata.vertices(:,1),Ctetdata.vertices(:,2),Ctetdata.vertices(:,3),100,colors,'filled');
end

% compute current based on voltages. GRAD(V)/R = I
bulkGradList = [ones(4*Ctetdata.numTetrahedra,1) reshape([Ctetdata.vertices(Ctetdata.tetrahedra(:,1),:) Ctetdata.vertices(Ctetdata.tetrahedra(:,2),:) Ctetdata.vertices(Ctetdata.tetrahedra(:,3),:) Ctetdata.vertices(Ctetdata.tetrahedra(:,4),:)]',3,[])'];
II = repmat([1:4*Ctetdata.numTetrahedra]',4,1);
Jt = reshape(repmat((0:(Ctetdata.numTetrahedra-1))*4,4,1),1,[]);
JJ = [Jt+1 Jt+2 Jt+3 Jt+4];
KK = bulkGradList(:);
bulkGradMat = sparse(II,JJ,KK,4*Ctetdata.numTetrahedra,4*Ctetdata.numTetrahedra);
bulkVoltages = reshape(phi(Ctetdata.tetrahedra)',[],1);
bulkAffineGrads = reshape(bulkGradMat\bulkVoltages,4,[])';
bulkLinearGrads = bulkAffineGrads(:,2:4);

if debugging
    scl = 10000;
    base = Ctetdata.tetBarycenters;
    head = bulkLinearGrads./norms(bulkLinearGrads,2,2);
    figure; hold all; axis equal; rotate3d on;
    quiver3(base(:,1),base(:,2),base(:,3),head(:,1),head(:,2),head(:,3))
    scatter3(Ctetdata.vertices(randomSource,1),Ctetdata.vertices(randomSource,2),Ctetdata.vertices(randomSource,3),500,'g.');
    scatter3(Ctetdata.vertices(randomSink,1),Ctetdata.vertices(randomSink,2),Ctetdata.vertices(randomSink,3),500,'r.');
end







