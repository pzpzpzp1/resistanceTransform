addpath('point2trimesh');
addpath('inpolyhedron');
clear all; close all;

%% declare misc parameters
resistivity = 1;
I = 1;
debugging = 0;
resolution = 20; % per edge
nMeasurements = 1000;
samplesizePerIter = 1000;
subdivide = false;

%% load random surface mesh
files = dir('../../10k_surface/');
rnum = randi(numel(files)-2)+2;
% rnum=7064; % standing ball thing.
rnum=4235;
filename = files(rnum).name;
[~,fname,ext] = fileparts(filename);

[V,F]=readOBJ(['../../10k_surface/' filename]);

%% build bounding box
margin = 1.2;
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

solutionVoltagesMat = nan(HMesh.nverts,nMeasurements);
Results = cell(nMeasurements,1);
inds = 2:HollowHMesh.nverts; % remove one variable to make laplacian full rank

bvp = boundaryVertPairs(measurements,:);
injectedCurrentMat = sparse(repmat(1:size(bvp,1),1,2), bvp(:), [ones(1,size(bvp,1)) -ones(1,size(bvp,1))], nMeasurements, HMesh.nverts)';

solutionVoltagesBulk = electricHollowLaplacian(inds,inds)\full(injectedCurrentMat(inds,:));
solutionVoltagesBulk = [zeros(1,size(solutionVoltagesBulk,2));solutionVoltagesBulk];

%% Move solved data onto full Box
injectedCurrentFull = sparse(HMesh.nverts, nMeasurements);
for measurementInd = 1:numel(measurements)
    measurement = measurements(measurementInd);
    sourcePos = boundaryVertPairs(measurement,1);
    sinkPos = boundaryVertPairs(measurement,2);
    
    % solve
    solutionVoltages = solutionVoltagesBulk(:,measurementInd);
    solutionVoltages = solutionVoltages-min(solutionVoltages); 
    
    if debugging && false
        sourceP = HollowHMesh.V2P(sourcePos,:);
        sinkP = HollowHMesh.V2P(sinkPos,:);
        rvs = HollowHMesh.V2P; colors = solutionVoltages/max(solutionVoltages);
        
        figure; hold all; rotate3d on; xlabel('voltages');
        scatter3(rvs(:,1),rvs(:,2),rvs(:,3),200,colors,'.');
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
    injectedCurrentFull(stind,measurementInd) = [-I,I];
    Results{measurementInd}.injectedCurrent = injectedCurrentFull(:,measurementInd);
        
    if debugging && false
        stinds = [Results{measurementInd}.sPosInd Results{measurementInd}.tPosInd];
        figure; hold all; rotate3d on; xlabel('voltages');
        Xs = HMesh.V2P(ridx,:);
        colors = measuredVoltages(ridx); colors = colors-min(colors); colors = colors/max(colors);
        scatter3(Xs(:,1),Xs(:,2),Xs(:,3),200,colors,'.');
        scatter3(HMesh.V2P(stinds,1),HMesh.V2P(stinds,2),HMesh.V2P(stinds,3),200,[0 1 0; 1 0 0],'filled');
    end
end




%% Solve for conductances using measured data!
% ground truth conductances
[ids, dists] = knnsearch(HollowHMesh.edgeCenters, HMesh.edgeCenters, 'K', 1);
conductances = (dists<1e-12)/resistivity;
conductancesGT = conductances; % ground truth conductance values on the full box. % one problem is that this makes the laplacian super degenerate.
electricLaplacianGT = HMesh.gradientOp'*diag(sparse(conductances))*HMesh.gradientOp;

% initialize variables
v0 = zeros(HMesh.nverts,nMeasurements);
conductances0 = conductancesGT*0+1;
conductances0 = conductances0/sum(conductances0)*sum(conductancesGT);
vol0 = sum(conductancesGT); % fixed volume

%% verify that ground truth conductances admits 0 energy solution
%{
    cvx_begin
        cvx_solver mosek
        variable v(HMesh.nverts,nMeasurements);
        
        % compute physical feasibility energy
        physFeas = norm(HMesh.gradientOp'*diag(sparse(conductancesGT))*HMesh.gradientOp ...
                *v-injectedCurrentFull);
        
        % compute proximity to measured values
        prox2Empirical = norm(v(bulkMeasuredVoltageIndices)-solutionVoltagesMat(bulkMeasuredVoltageIndices));
        minimize prox2Empirical + alpha*physFeas;
    cvx_end
%}

%% minimize biconvex energy by alternating
figure; hold all; rotate3d on; xlabel('GT conductance');
ptc2 = patch('Faces',HollowHMesh.F2V(HollowHMesh.isBoundaryFace,:),'Vertices',HollowHMesh.V2P,'FaceColor','green','EdgeColor','none'); alpha(ptc2,.1)
xlim(BB(:,1)'+[1 -1]*1e-1);ylim(BB(:,2)'+[1 -1]*1e-1);zlim(BB(:,3)'+[1 -1]*1e-1);
thresh = find(1-conductancesGT > 1e-2);
scatter3(HMesh.edgeCenters(thresh,1),HMesh.edgeCenters(thresh,2),HMesh.edgeCenters(thresh,3),5,conductancesGT(thresh));
converged = false;
alphac = 1;
f1 = figure; hold all; rotate3d on; sctr = scatter3([],[],[]); xlabel('Iterated Conductances');
ptc2 = patch('Faces',HollowHMesh.F2V(HollowHMesh.isBoundaryFace,:),'Vertices',HollowHMesh.V2P,'FaceColor','green','EdgeColor','none'); alpha(ptc2,.1)
xlim(BB(:,1)'+[1 -1]*1e-1);ylim(BB(:,2)'+[1 -1]*1e-1);zlim(BB(:,3)'+[1 -1]*1e-1);
counter = 1;
while ~converged
    fprintf('Iteration %d\n',counter);
    selectedMeasurements = randsample(nMeasurements,samplesizePerIter,false);
    
    %% alternate minimizing |Lv-I|^2 + alpha*|v_m-v_e|^2
    t1=tic;
    
    % fix conductances, minimize w.r.t. voltages
    vmeasured = solutionVoltagesMat(HMesh.isBoundaryVerts, selectedMeasurements);
    Imeasured = injectedCurrentFull(:,selectedMeasurements);
    Lap = HMesh.gradientOp'*diag(sparse(conductances0))*HMesh.gradientOp;
    VM = sparse(1:sum(HMesh.isBoundaryVerts),find(HMesh.isBoundaryVerts),ones(sum(HMesh.isBoundaryVerts),1),sum(HMesh.isBoundaryVerts),numel(HMesh.isBoundaryVerts)); %virtual measurement operator
    v = (Lap'*Lap + alphac * VM'*VM)\(Imeasured' * Lap + alphac * vmeasured'*VM)';
    
    % fix voltages, minimize relative to conductances
    cvx_begin
        cvx_solver mosek
        variable conductances0(HMesh.nedges,1);
        
        % regularize conductance?
        
        % compute physical feasibility energy
        physFeas = norm(HMesh.gradientOp' * diag(sparse(conductances0)) * HMesh.gradientOp ...
                *v - injectedCurrentFull(:,selectedMeasurements),'fro');
            
        minimize physFeas;
        subject to
            conductances0 <= 1;
            conductances0 >= 0;
            sum(conductances0) == vol0;
            conductances0(HMesh.isBoundaryEdge)==ones(sum(HMesh.isBoundaryEdge),1)
    cvx_end
    telapsed(counter) = toc(t1);
    conductancesPrev = conductances0;
    
    if debugging
        figure(f1); 
        delete(sctr);
        thresh = find(1-conductancesPrev > 1e-2);
        sctr = scatter3(HMesh.edgeCenters(thresh,1),HMesh.edgeCenters(thresh,2),HMesh.edgeCenters(thresh,3),5,conductancesPrev(thresh));
        xlabel(sprintf('Iterated Conductances %d',counter));
    end
    deviation(counter) = norm(conductancesPrev-conductancesGT);
    counter = counter + 1;    
end

%% threshold based on volume constraint.
[sortedC, sortedInds] = sort(conductancesPrev);
keepinds = sortedInds(1:(numel(conductances0)-vol0));
figure; hold all; rotate3d on; 
xlim(BB(:,1)'+[1 -1]*1e-1);ylim(BB(:,2)'+[1 -1]*1e-1);zlim(BB(:,3)'+[1 -1]*1e-1);
scatter3(HMesh.edgeCenters(keepinds,1),HMesh.edgeCenters(keepinds,2),HMesh.edgeCenters(keepinds,3),30,ones(numel(keepinds),1),'green','filled');
thresholdedConductances = ones(numel(conductances0),1); thresholdedConductances(keepinds)=0; showinds = find((thresholdedConductances-conductancesGT)~=0)
scatter3(HMesh.edgeCenters(showinds,1),HMesh.edgeCenters(showinds,2),HMesh.edgeCenters(showinds,3),5,thresholdedConductances(showinds)-conductancesGT(showinds));
title('Diff with GT'); xlabel('Green is the resulting reconstruction. Yellow is the missing parts. Blue are the extra parts.');

% figure; semilogy(deviation); %ylim([0 max(deviation)]);
% figure; plot(telapsed);

















