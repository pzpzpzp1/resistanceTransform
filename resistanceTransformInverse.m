addpath('point2trimesh');

%% declare misc parameters
resistivity = 1;
I = 1;
debugging = 1;
premade = 1;
resolution = 20; % 30 per edge
nMeasurements = 500;
subdivide = false;

%% load tet mesh
if premade
    files = dir('meshes/vtkComplementTetMeshes');
    if numel(files)==2
        premade = false;
    else
        filename = files(randi(numel(files)-2)+2).name;
        [ext,fnamec,ext] = fileparts(filename);
    end
end
if ~premade
    files = dir('../../10k_tetmesh/');
    filename = files(randi(numel(files)-2)+2).name;
    [ext,fname,ext] = fileparts(filename);
    fnamec = [fname '_c'];

    generateComplementMesh(fname);
end

%% LOAD COMPLEMENT TET MESH. Possibly with subdivision.
fid = fopen(['meshes/vtkComplementTetMeshes/' fnamec '.vtk']);
[CVerts, CTets] = loadVTKTET(fid);
fclose(fid);

if subdivide
    [CVerts,CTets] = subdivideTetMesh(CVerts,CTets);
end

tetra = triangulation(CTets, CVerts);
[TetL, TetM] = GeometricPrimalLM(tetra); % laplacian and mass matrix
Ctetdata = getTetDataRT(CTets,CVerts,1);

%% Extract rectangle mesh and rectangle mesh properties
BB = [min(CVerts); max(CVerts)];
xverts = linspace(BB(1,1),BB(2,1),resolution)';
yverts = linspace(BB(1,2),BB(2,2),resolution);
zverts = permute(linspace(BB(1,3),BB(2,3),resolution),[1 3 2]);
X = repmat(xverts,1,resolution,resolution);
Y = repmat(yverts,resolution,1,resolution);
Z = repmat(zverts,resolution,resolution,1);

indices = reshape(1:numel(X),[resolution resolution resolution]);
rectVerts = [X(:) Y(:) Z(:)];
clear rectCells;
rectCells(:,:,:,1) = indices((1:resolution-1),(1:resolution-1),(1:resolution-1));
rectCells(:,:,:,2) = indices((1:resolution-1),(1:resolution-1),(1:resolution-1)+1);
rectCells(:,:,:,3) = indices((1:resolution-1),(1:resolution-1)+1,(1:resolution-1));
rectCells(:,:,:,4) = indices((1:resolution-1),(1:resolution-1)+1,(1:resolution-1)+1);
rectCells(:,:,:,5) = indices((1:resolution-1)+1,(1:resolution-1),(1:resolution-1));
rectCells(:,:,:,6) = indices((1:resolution-1)+1,(1:resolution-1),(1:resolution-1)+1);
rectCells(:,:,:,7) = indices((1:resolution-1)+1,(1:resolution-1)+1,(1:resolution-1));
rectCells(:,:,:,8) = indices((1:resolution-1)+1,(1:resolution-1)+1,(1:resolution-1)+1);
rectCells = reshape(permute(rectCells,[4 1 2 3]),8,[])';

if debugging && false
    figure; hold all; axis equal; rotate3d on;
    scatter3(rectVerts(:,1),rectVerts(:,2),rectVerts(:,3),'.')
    
    figure; hold all; axis equal; rotate3d on;
    for i=1:size(rectCells,1)
        clc; hold all; axis equal; rotate3d on;
        verts = rectVerts(rectCells(i,:),:);
        scatter3(verts(:,1),verts(:,2),verts(:,3))
        pause;
    end
end

rectEdges = unique(sort(reshape(rectCells(:,[1 2 2 4 4 3 3 1 [1 2 2 4 4 3 3 1]+4 1 5 2 6 3 7 4 8])',2,[])',2),'rows');
nE = size(rectEdges,1);
nV = size(rectVerts,1);
nC = size(rectCells,1);
rectIncidenceMatrix = sparse([1:nE 1:nE]', [rectEdges(:,1);rectEdges(:,2)], [ones(nE,1) -ones(nE,1)], nE, nV);
edgelengths = norms(rectVerts(rectEdges(:,1),:)-rectVerts(rectEdges(:,2),:),2,2);
gradientOp = rectIncidenceMatrix ./ edgelengths;
vertDegree = full(sum(abs(rectIncidenceMatrix)))';
BoundaryVertices = find(vertDegree~=6);
nBv = numel(BoundaryVertices);
II = repmat(BoundaryVertices,nBv,1);
JJ = repmat(BoundaryVertices',nBv,1); JJ = JJ(:);
boundaryVertPairs = [II JJ];
nBvPairs = size(boundaryVertPairs,1);

%% In silico data collection: 
% some measurements can fail if pointlocation floating precision isn't good.
% successful measurements can contain nans, for the same reason.
measurements = randsample(nBvPairs,nMeasurements,false);

Results = cell(numel(measurements),1);
exceptions = {};
successfulMeasurements = zeros(nMeasurements,1);
for measurementInd = 1:numel(measurements)
    measurement = measurements(measurementInd);
    sourcePos = rectVerts(boundaryVertPairs(measurement,1),:);
    sinkPos = rectVerts(boundaryVertPairs(measurement,2),:);

    try
        measuredVoltages = measureVoltages(Ctetdata, tetra, TetL, sourcePos, sinkPos, I, rectVerts(BoundaryVertices,:), resistivity);
    catch exception
        exceptions{end+1} = exception;
        continue;
    end

    Results{measurementInd}.measuredVoltages = measuredVoltages;
    Results{measurementInd}.sPosInd = boundaryVertPairs(measurement,1);
    Results{measurementInd}.tPosInd = boundaryVertPairs(measurement,2);
    successfulMeasurements(measurementInd)=1;
end
nSuccess = nMeasurements-numel(exceptions);

%% Solve for conductances using measured data!





















