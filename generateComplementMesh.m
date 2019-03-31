function generateComplementMesh(fname)
    debugging = 0;

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
        len = norm(BB(2,:)-BB(1,:))*.05; % default of tetwild
        len = min(BB(2,:)-BB(1,:))/20;
        system(sprintf('TetWild.exe %s -q --is-quiet --output %s.msh -l %f',surfaceFname,gmshTetMesh,len));
        delete meshes/gmshTetMeshes/*_sf.obj;
        
        % gmsh
        system(sprintf('%s %s.msh -save -o %s.vtk',gmshbase, gmshTetMesh,['meshes/vtkComplementTetMeshes/' fnamec]));
        
    end
end