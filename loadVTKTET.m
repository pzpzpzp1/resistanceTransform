% returns vertex list and tets 1-indexed into verts.
function [X,T] = loadVTKTET(fid)

    line = ''; safety = 1;
    while ~contains(line,'POINTS') && safety < 10
        safety = safety + 1;
        line=fgets(fid);
    end

%     textscan(fid, '# vtk DataFile Version 2.0');
%     textscan(fid, 'TET');
%     textscan(fid, 'tetra');
%     textscan(fid, 'ASCII');
%     textscan(fid, '');
%     textscan(fid, 'DATASET UNSTRUCTURED_GRID');
% 
%     NV = textscan(fid, 'POINTS %d float');
%     NV=NV{1};

    V = textscan(fid, '%f %f %f', 'CollectOutput', true, 'MultipleDelimsAsOne', true);
    X=V{1};

    textscan(fid, 'CELLS %d %d');
    TC = textscan(fid, '4  %d %d %d %d', 'CollectOutput', true, 'MultipleDelimsAsOne', true);
    T=TC{1}; % 0 indexed

    T=double(T+1);
end