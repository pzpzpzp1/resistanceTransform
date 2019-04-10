addpath('point2trimesh');
addpath('inpolyhedron');
clear all; close all;

%% declare misc parameters
resistivity = 1;
I = 1;
debugging = 1;
resolution = 30; % per edge
nMeasurements = 500;
subdivide = false;

%% load random surface mesh
files = dir('../../10k_surface/');
rnum = randi(numel(files)-2)+2;
filename = files(rnum).name;
[~,fname,ext] = fileparts(filename);

[V,F]=readOBJ(['../../10k_surface/' filename]);

%% build bounding box
margin = 1.1;
means = (max(V)+min(V))/2;
BB = [means-(means-min(V))*margin; (max(V)-means)*margin+means];
assert(all(BB(1,:)<min(V) & BB(2,:)>max(V)));

%% Extract base rectangle mesh: (rectVerts, rectCells)
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
rectCells = rectCells(:,[1 2 4 3 5 6 8 7]); % permute into standard hmesh form
HMesh = loadHexMesh(rectVerts, rectCells); % full box.

%% Find which box barycenters are not in the original Obj file. Construct conductance matrix per edge of original mesh.
rectBarycenters = (rectVerts(rectCells(:,1),:)+...
                    rectVerts(rectCells(:,2),:)+...
                    rectVerts(rectCells(:,3),:)+...
                    rectVerts(rectCells(:,4),:)+...
                    rectVerts(rectCells(:,5),:)+...
                    rectVerts(rectCells(:,6),:)+...
                    rectVerts(rectCells(:,7),:)+...
                    rectVerts(rectCells(:,8),:))/8;
cellsToRemove = inpolyhedron(F,V,rectBarycenters);

rectCells(cellsToRemove,:)=[];
[rectVertsO, rectCellsO, idx] = minimizeHexMesh(rectVerts, rectCells); 
HollowHMesh = loadHexMesh(rectVertsO, rectCellsO);

if debugging
    figure; 
    subplot(1,2,1); hold all; axis equal; rotate3d on; xlabel('points in object');
    %scatter3(rectBarycenters(cellsToRemove,1),rectBarycenters(cellsToRemove,2),rectBarycenters(cellsToRemove,3),10,'r','filled');
    ptc1 = patch('Faces',F,'Vertices',V,'FaceColor','green','EdgeColor','none'); alpha(ptc1,.1);
    
    subplot(1,2,2); hold all; axis equal; rotate3d on; xlabel('voxel representation');
    ptc2 = patch('Faces',HollowHMesh.F2V(HollowHMesh.isBoundaryFace,:),'Vertices',HollowHMesh.V2P,'FaceColor','green','EdgeColor','black');
    xlim(BB(:,1)'+[1 -1]*1e-1);
    ylim(BB(:,2)'+[1 -1]*1e-1);
    zlim(BB(:,3)'+[1 -1]*1e-1);
    alpha(ptc2,.9)
end

%% Get Data from forward simulation
HollowHMesh.isBoundingBoxVerts = HMesh.isBoundaryVerts(idx);
electricHollowLaplacian = HollowHMesh.gradientOp'*HollowHMesh.gradientOp;
[ridx, dists] = knnsearch(HMesh.V2P, HollowHMesh.V2P(HollowHMesh.isBoundingBoxVerts,:), 'K', 1);
assert(norm(dists)==0);

BoundaryVertices = find(HollowHMesh.isBoundingBoxVerts);
nBv = numel(BoundaryVertices);
II = repmat(BoundaryVertices,nBv,1);
JJ = repmat(BoundaryVertices',nBv,1); JJ = JJ(:);
boundaryVertPairs = [II JJ];
boundaryVertPairs = boundaryVertPairs(boundaryVertPairs(:,1)~=boundaryVertPairs(:,2),:);
nBvPairs = size(boundaryVertPairs,1);
assert(nMeasurements <= nBvPairs); % can't take more samples than exist with a certain discretization
measurements = randsample(nBvPairs,nMeasurements,false);

injectedCurrentMat = sparse(HMesh.nverts,nMeasurements);
solutionVoltagesMat = nan(HMesh.nverts,nMeasurements);
Results = cell(nMeasurements,1);
inds = 2:HollowHMesh.nverts; % remove one variable to make laplacian full rank
L=chol(electricHollowLaplacian(inds,inds)); % chol factor to speed up solve. ~2x speedup per solve, for which there are many!
for measurementInd = 1:numel(measurements)
    fprintf('Made measurement %d out of %d\n', measurementInd, nMeasurements);
    measurement = measurements(measurementInd);
    sourcePos = boundaryVertPairs(measurement,1);
    sinkPos = boundaryVertPairs(measurement,2);
    
    % construct current injection vector. Technically integrating a scaled delta function with the hex basis.
    currentVector = zeros(HollowHMesh.nverts,1);
    currentVector(sinkPos)=-I;
    currentVector(sourcePos)=I;
    
    % solve
    % solutionVoltages = [0; electricHollowLaplacian(inds,inds)\currentVector(inds)]; % This is the expensive version that isn't used.
    solutionVoltages = [0; L\(L'\currentVector(inds))];
    solutionVoltages = solutionVoltages-min(solutionVoltages); 
    
    if debugging && false
        sourceP = HollowHMesh.V2P(sourcePos,:);
        sinkP = HollowHMesh.V2P(sinkPos,:);
        rvs = HollowHMesh.V2P; colors = solutionVoltages/max(solutionVoltages);
        
        figure; hold all; rotate3d on; xlabel('voltages');
        scatter3(rvs(:,1),rvs(:,2),rvs(:,3),200,[colors colors colors],'.');
        scatter3(sourceP(:,1),sourceP(:,2),sourceP(:,3),300,'r','filled');
        scatter3(sinkP(:,1),sinkP(:,2),sinkP(:,3),300,'g','filled');
        xlim(BB(:,1)');
        ylim(BB(:,2)');
        zlim(BB(:,3)');
    end
    
    % measure at boundary: only the boundary part of the solution can be 'measured'
    measuredVoltages = nan(HMesh.nverts,1);
    measuredVoltages(ridx) = solutionVoltages(HollowHMesh.isBoundingBoxVerts);
    Results{measurementInd}.sPosInd = ridx(find(sinkPos==find(HollowHMesh.isBoundingBoxVerts)));
    Results{measurementInd}.tPosInd = ridx(find(sourcePos==find(HollowHMesh.isBoundingBoxVerts)));
    Results{measurementInd}.measuredVoltages = measuredVoltages; 
    solutionVoltagesMat(ridx,measurementInd) = solutionVoltages(HollowHMesh.isBoundingBoxVerts);
    stind = [Results{measurementInd}.sPosInd Results{measurementInd}.tPosInd];
    injectedCurrentMat(stind,measurementInd) = [-I,I];
    Results{measurementInd}.injectedCurrent = injectedCurrentMat(:,measurementInd);
    
    if debugging && false
        stinds = [Results{measurementInd}.sPosInd Results{measurementInd}.tPosInd];
        figure; hold all; rotate3d on; xlabel('voltages');
        Xs = HMesh.V2P(ridx,:);
        colors = measuredVoltages(ridx); colors = colors-min(colors); colors = colors/max(colors);
        scatter3(Xs(:,1),Xs(:,2),Xs(:,3),200,[colors colors colors],'.');
        scatter3(HMesh.V2P(stinds,1),HMesh.V2P(stinds,2),HMesh.V2P(stinds,3),200,[0 1 0; 1 0 0],'filled');
    end
end

%% Solve for conductances using measured data!
% ground truth conductances
[ids, dists] = knnsearch(HollowHMesh.edgeCenters, HMesh.edgeCenters, 'K', 1);
conductances = (dists<1e-12)/resistivity;
conductancesGT = conductances; % ground truth conductance values on the full box. % one problem is that this makes the laplacian super degenerate.
electricLaplacianGT = HMesh.gradientOp'*diag(sparse(conductances))*HMesh.gradientOp;
bulkMeasuredVoltageIndices = find(~isnan(solutionVoltagesMat));


% initialize variables
v0 = zeros(HMesh.nverts,nMeasurements);
conductances0 = conductancesGT*0+1;

converged = false;
while ~converged
    %% solve for v0 holding constant conductances
    % minimize voltage difference from empirical
    % subject to current injection constraints
    cvx_begin
        cvx_solver mosek
        
    cvx_end
    
    
    
    
    %% solve for conductances, holding voltages
    % minimize total variation of conductance
    % subject to current injection constraints
    % subject to volume constraint
    cvx_begin
        cvx_solver mosek
        
    cvx_end
    
    
    
end

% quantity to minimize
measuredVoltages = v0*nan;
measuredVoltages(BoundaryVertices,:)=measuredBoundaryVoltages;
measureInds = find(~isnan(measuredVoltages));
flattenedVoltages = measuredVoltages(measureInds);
% v0(measureInds) == flattenedVoltages;

% constraints
% for i = 1:nSuccess
%     rect.gradient' * spdiag(conductances0) * rect.gradient * v0(:,i) == currentInjections(:,i);
% end




















