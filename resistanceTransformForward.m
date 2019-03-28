

%% declare misc parameters
resistivity = 1;
va = 0;
vb = 1;


%% load tet mesh
files = dir('../../10k_tetmesh');
filename = files(randi(numel(files)-2)+2).name;
[ext,fname,ext] = fileparts(filename);

if exist(['meshes/' fname '.vtk'])~=2
    if ispc
        [ret, name] = system('hostname');
        if contains(name,'DESKTOP-UH1G0VB')
            system(sprintf('D:/Documents/gmsh-4.2.2-Windows64/gmsh.exe %s.msh -save -o %s.vtk',['../../10k_tetmesh/' fname],['meshes/' fname]));
        else
            error('new computer?');
        end
    else
        error('new computer?');
    end
end
fid = fopen(['meshes/' fname '.vtk']);
[X,T] = loadVTKTET(fid);
fclose(fid);

%% build bounding box

means = [max(X)+min(X)]/2;
BB = [means-(means-min(X))*1.1; (max(X)-means)*1.1+means];
assert(all(BB(1,:)<min(X) & BB(2,:)>max(X)));

%% create complement tet mesh

%% load complement tet mesh

%% Solve poisson problem for electrical properties




