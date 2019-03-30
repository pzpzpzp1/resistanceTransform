
%% declare misc parameters
resistivity = 1;
I = 1;
debugging = 1;

%% load tet mesh
files = dir('../../10k_tetmesh');
filename = files(randi(numel(files)-2)+2).name;
[ext,fname,ext] = fileparts(filename);
fnamec = [fname '_c'];

if exist(['meshes/vtkTetMeshes/' fnamec '.vtk'])~=2
    if ispc
        [ret, name] = system('hostname');
        if contains(name,'DESKTOP-UH1G0VB')
            gmshbase = 'D:/Documents/gmsh-4.2.2-Windows64/gmsh.exe';
        elseif contains(name,'DESKTOP-5UNMMEJ')
            gmshbase = 'C:\Users\pzpzp\Documents\gmsh-4.2.2-Windows64\gmsh.exe';
        else
            error('new computer?');
        end
    else
        error('new computer?');
    end
    system(sprintf('%s %s.msh -save -o %s.vtk',gmshbase,['../../10k_tetmesh/' fname],['meshes/vtkTetMeshes/' fname]));
            
    
    fid = fopen(['meshes/vtkTetMeshes/' fname '.vtk']);
    [Verts,Tets] = loadVTKTET(fid);
    fclose(fid);
    if size(Verts,1)==0
        throw(['Empty Mesh ' filename]);
    end

    %% build bounding box
    margin = 1.3;
    means = (max(Verts)+min(Verts))/2;
    BB = [means-(means-min(Verts))*margin; (max(Verts)-means)*margin+means];
    assert(all(BB(1,:)<min(Verts) & BB(2,:)>max(Verts)));

    %% create complement surface mesh
    % extract surface triangle mesh. throw if nonmanifold
    tetdata = getTetDataRT(Tets,Verts,1);
    BoundaryTriangles = tetdata.triangles(tetdata.isBoundaryTriangle==1,:);
    [SurfaceVerts, SurfaceTriangles]=minimizeMesh(Verts,BoundaryTriangles);
    SurfaceTriangles = SurfaceTriangles(:,[1 2 3]); % flip or unflip boundary inside out. seems unnecessary actually...

    % add on bounding box mesh
    rectIndices = dec2bin(0:7)-47;
    rectVerts = [BB(rectIndices(:,1),1) BB(rectIndices(:,2),2) BB(rectIndices(:,3),3)];
    rectK = convhull(rectVerts(:,1),rectVerts(:,2),rectVerts(:,3));

    % combine surfaces
    newSurfaceX = [SurfaceVerts; rectVerts];
    newSurfaceTriangles = [SurfaceTriangles; rectK+size(SurfaceVerts,1)];

    if debugging
        figure; hold all; axis equal; rotate3d on;
        ptc=patch('Faces',newSurfaceTriangles,'Vertices',newSurfaceX,'FaceColor','green','EdgeColor','black');
        alpha(ptc,.1);
    end

    %% create complement tet mesh
    % save complement surface
    surfaceFname = ['meshes/surfaceOffMeshes/'  fname '.off'];
    gmshTetMesh = ['meshes/gmshTetMeshes/'  fnamec];
    writeOff(surfaceFname, newSurfaceX, newSurfaceTriangles);

    % call tetwild
    system(sprintf('TetWild.exe %s -q --is-quiet --output %s.msh',surfaceFname,gmshTetMesh));
    

    % gmsh
    system(sprintf('%s %s.msh -save -o %s.vtk',gmshbase, gmshTetMesh,['meshes/vtkTetMeshes/' fnamec]));
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







