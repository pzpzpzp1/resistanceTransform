
%% declare misc parameters
resistivity = 1;
va = 0;
vb = 1;
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
        else
            error('new computer?');
        end
    else
        error('new computer?');
    end

    fid = fopen(['meshes/vtkTetMeshes/' fname '.vtk']);
    [Verts,Tets] = loadVTKTET(fid);
    fclose(fid);

    %% build bounding box
    means = (max(Verts)+min(Verts))/2;
    BB = [means-(means-min(Verts))*1.1; (max(Verts)-means)*1.1+means];
    assert(all(BB(1,:)<min(Verts) & BB(2,:)>max(Verts)));

    %% create complement surface mesh
    % extract surface triangle mesh. throw if nonmanifold
    tetdata = getTetDataRT(Tets,Verts,1);
    BoundaryTriangles = tetdata.triangles(tetdata.isBoundaryTriangle==1,:);
    [SurfaceVerts, SurfaceTriangles]=minimizeMesh(Verts,BoundaryTriangles);
    SurfaceTriangles = SurfaceTriangles(:,[1 3 2]); % flip boundary inside out.

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
    system(sprintf('TetWild.exe %s -q --is-quiet --output %s',surfaceFname,gmshTetMesh));

    % gmsh
    system(sprintf('%s %s.msh -save -o %s.vtk',gmshbase, [gmshTetMesh '.msh'],['meshes/vtkTetMeshes/' fnamec]));
end

%% LOAD COMPLEMENT TET MESH
fid = fopen(['meshes/vtkTetMeshes/' fnamec '.vtk']);
[CVerts, CTets] = loadVTKTET(fid);
fclose(fid);

tetra = triangulation(CTets, CVerts);
[TetL, TetM] = GeometricPrimalLM(tetra); % laplacian and mass matrix

%% Solve poisson problem for electrical properties




