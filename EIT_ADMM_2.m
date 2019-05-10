addpath('point2trimesh');
addpath('inpolyhedron');
clear all; close all;
profile on;

%% declare misc parameters
verifySubparts = 1;
resistivity = 1;
I = 1;
debugging = 0;
resolution = 20; % per edge
nMeasurements = 500;
samplesizePerIter = nMeasurements;
subdivide = false;

%% load random surface mesh
files = dir('../../10k_surface/');
rnum = randi(numel(files)-2)+2;
% rnum=7064; % standing ball thing.
% rnum=4235;
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
volt0 = zeros(HMesh.nverts,nMeasurements);
conductances0 = conductancesGT*0+1;
conductances0 = conductances0/sum(conductances0)*sum(conductancesGT);
v0 = sum(conductancesGT); % fixed volume




%% ADMM section
% initialize
boundaryEdgeInds = find(HMesh.isBoundaryEdge);
interiorEdgeInds = find(HMesh.isBoundaryEdge);
boundaryInds = find(HMesh.isBoundaryVerts);
interiorInds = find(~HMesh.isBoundaryVerts);
M = sparse(1:numel(boundaryInds),boundaryInds,boundaryInds*0+1,numel(boundaryInds),HMesh.nverts);
D = HMesh.gradientOp;
phi = volt0;
phihat = solutionVoltagesMat(boundaryInds, :);
sigmavec = conductances0 + rand(size(conductances0));
sigma = sparse(diag(sigmavec));
sigmaT = sigma + diag(rand(size(sigma,1),1)*2);
J = injectedCurrentFull;
Z = randn(size(J))*10;
w = randn*10;
z = rand(size(sigma,1),1)-.5;
rho1 = rand*100;
rho2 = rand*100;

%% verify that ground truth conductances admits 0 energy solution
sigmaGT = sparse(diag(conductancesGT));
phiGT = (D'*sigmaGT*D\J);
E_GT = norm(full(D'*sigmaGT*D*phiGT-J));
assert(E_GT <= 1e-8);

%% solve ADMM iterations
converged = false;
counter = 1;
energies = [];
times = [];
while ~converged
    fprintf('Iteration %d\n',counter);
    % augmented lagrangian energy
    energies(counter) = norm(phihat-M*phi,'fro')^2 + trace((D'*sigma*D*phi-J)*Z') + w*(sum(diag(sigma))-v0)...
        + rho1/2*norm(D'*sigma*D*phi-J,'fro')^2 + rho2/2*(sum(diag(sigma))-v0)^2 + z'*diag(sigmaT-sigma) + rho1/2*norm(diag(sigmaT-sigma))^2;
    tperiter = tic;
    
    % dual update
    Z = Z + rho1*(D'*sigma*D*phi-J);
    w = w + rho2*(sum(sum(sigma))-v0);
    z = z + rho1*diag(sigmaT-sigma);
    
    % update phi
    L = D'*sigma*D;
    S = 2*(M'*M) + rho1*(L'*L);
    t = (2*phihat'*M - Z'*L + rho1*J'*L)';
    phi = S\t;
    
    % update sigmaT
    sigmaT = sigma-diag(z)/rho1;
    sigmaT(sigmaT>1) = 1;
    sigmaT(sigmaT<0) = 0;
    
    if verifySubparts && false
        % sigmaT
        cvx_begin
            cvx_precision best;
            cvx_solver mosek;
            variable sigmaTm(size(sigmaT));
            z1 = z'*diag(sigmaTm-sigma) + rho1/2*pow_pos(norm(diag(sigmaTm-sigma),'fro'),2);
            minimize z1
            subject to 
                diag(sigmaTm) <= 1
                diag(sigmaTm) >= 0
                sigmaTm == diag(diag(sigmaTm));
        cvx_end
        z2 = z'*diag(sigmaT-sigma) + rho1/2*pow_pos(norm(diag(sigmaT-sigma),'fro'),2);
        assert(max(max(abs(z1-z2)))<1e-2);

        % phi
        cvx_begin
            cvx_precision best;
            cvx_solver mosek;
            variable phim(size(phi));
            z1 = pow_pos(norm(phihat-M*phim,'fro'),2) + trace((D'*sigma*D*phim-J)*Z') + w*(sum(diag(sigma))-v0)...
                    + rho1/2*pow_pos(norm(D'*sigma*D*phim-J,'fro'),2) + rho2/2*pow_abs(sum(diag(sigma))-v0,2);
            minimize z1
        cvx_end
        z2 = pow_pos(norm(phihat-M*phi,'fro'),2) + trace((D'*sigma*D*phi-J)*Z') + w*(sum(diag(sigma))-v0)...
        + rho1/2*pow_pos(norm(D'*sigma*D*phi-J,'fro'),2) + rho2/2*pow_abs(sum(diag(sigma))-v0,2);
        assert(norm(z1-z2)<1e-3);
    end
    
    % update sigma
    Dt = D';
    Dphi = D*phi;
    %G = full(Dt).*permute(Dphi,[3 1 2]); % can't compute this. too dense and large
    %Kt = sum(sum(permute(G,[2 4 1 3]).*permute(G,[4 2 1 3]),3),4);
    %ht = sum(sum(G.*permute(full(J),[1 3 2]),1),3)';
    K = (Dphi*Dphi').*(Dt'*Dt);
    h = sum((Dt'*J).*Dphi,2);
    g = diag(D*Z*phi'*D') + w;
    A = rho1*K + rho2;
    b = (rho1*h' + rho2*v0*ones(size(h')) - g')';
    Aprime = A + rho1*speye(size(A));
    bprime = b + z + rho1*diag(sigmaT);
    %sigmavec = Aprime\bprime;
     knownInds = find(HMesh.isBoundaryEdge);
     unknownInds = find(~HMesh.isBoundaryEdge);
     sigmavec(knownInds,:) = 1;
     sigmavec(unknownInds,:) = Aprime(unknownInds,unknownInds)\(bprime(unknownInds,:) - sum(Aprime(knownInds,unknownInds))');
    sigma = diag(sparse(sigmavec));
    
    if verifySubparts && false
        %{
        G = full(Dt).*permute(Dphi,[3 1 2]);
        Efull0 = trace((D'*sigma*D*phi-J)*Z') + w*(sum(diag(sigma))-v0) + rho1/2*norm(D'*sigma*D*phi-J,'fro')^2 + rho2/2*(sum(diag(sigma))-v0)^2;
        Efull1 = trace((D*Z*phi'*D'+w)*sigma') + rho1/2*norm(D'*sigma*D*phi-J,'fro')^2 + rho2/2*(sum(diag(sigma))-v0)^2 + trace((-J)*Z') + w*(-v0);
        Efull2 = trace((D*Z*phi'*D'+w)*sigma') + rho1/2*norm(permute(sum(G.*full(diag(sigma))',2),[1 3 2])-J,'fro')^2 + rho2/2*(sum(diag(sigma))-v0)^2;
        Efull3 = trace((D*Z*phi'*D'+w)*sigma') + rho1/2*norm(permute(sum(G.*full(diag(sigma))',2),[1 3 2])-J,'fro')^2 + rho2/2*(ones(size(sigma,1),1)'*diag(sigma)-v0)^2;
        Efull4 = rho1/2*norm(permute(sum(G.*full(diag(sigma))',2),[1 3 2])-J,'fro')^2 + rho2/2*(ones(size(sigma,1),1)'*diag(sigma)-v0)^2;
        Gsigma = permute(sum(G.*full(diag(sigma))',2),[1 3 2]);
        Efull5 = rho1/2*norm(Gsigma-J,'fro')^2 + rho2/2*(ones(size(sigma,1),1)'*diag(sigma)-v0)^2;
        E51 = rho1/2*norm(Gsigma-J,'fro')^2;
        Efull6 = rho1/2*sum(sum(Gsigma.*Gsigma))-rho1*sum(sum(Gsigma.*J))+rho1/2*sum(sum(J.^2)) + rho2/2*(ones(size(sigma,1),1)'*diag(sigma)-v0)^2;
        Efull7 = rho1/2*sum(sum(Gsigma.*Gsigma))-rho1*sum(sum(Gsigma.*J)) + rho2/2*(ones(size(sigma,1),1)'*diag(sigma)-v0)^2;
        Efull8 = rho1/2*sum(sum(Gsigma.*Gsigma))-rho1*sum(sum(Gsigma.*J)) + rho2/2*diag(sigma)'*ones(size(sigma,1))*diag(sigma) - rho2*v0*sum(diag(sigma)) + rho2/2*v0^2;
        Efull9 = rho1/2*sum(sum(Gsigma.*Gsigma)) - rho1*sum(sum(Gsigma.*J)) + rho2/2*diag(sigma)'*ones(size(sigma,1))*diag(sigma) - rho2*v0*sum(diag(sigma));
        E91 = rho1*h'*sigmavec
        E92 = rho1*sum(sum(Gsigma.*J))
        Efull10 = rho1/2*sum(sum(Gsigma.*Gsigma)) - rho1*h'*sigmavec + rho2/2*diag(sigma)'*ones(size(sigma,1))*diag(sigma) - rho2*v0*sum(diag(sigma));
        E101 = rho1/2*sum(sum(Gsigma.*Gsigma));
        E102 = rho1/2*sigmavec'*K*sigmavec;
        Efull11 = rho1/2*sigmavec'*K*sigmavec - rho1*h'*sigmavec + rho2/2*diag(sigma)'*ones(size(sigma,1))*diag(sigma) - rho2*v0*sum(diag(sigma));
        Efull12 = 1/2*sigmavec'*A*sigmavec - b'*sigmavec;
        
        cvx_begin
            cvx_precision best;
            cvx_solver mosek;
            variable sigmam(size(sigma)) diagonal;
            sigmavecm = diag(sigmam);
            z = pow_pos(norm(phihat-M*phi,'fro'),2) + trace((D'*sigmam*D*phi-J)*Z') + w*(sum(diag(sigmam))-v0)...
                    + rho1/2*pow_pos(norm(D'*sigmam*D*phi-J,'fro'),2) + rho2/2*pow_abs(sum(diag(sigmam))-v0,2);
            minimize z            
        cvx_end
        sigmauncon = diag(A\b);
        z2 = pow_pos(norm(phihat-M*phi,'fro'),2) + trace((D'*sigmauncon*D*phi-J)*Z') + w*(sum(diag(sigmauncon))-v0)...
            + rho1/2*pow_pos(norm(D'*sigmauncon*D*phi-J,'fro'),2) + rho2/2*pow_abs(sum(diag(sigmauncon))-v0,2);
        
        cvx_begin
            cvx_precision best;
            cvx_solver mosek;
            variable xm(size(sigmavec));
            z3 = pow_pos(norm(C*xm-d,'fro'),2);
            minimize z3
            subject to
                xm >= 0
                xm <= 1
                sum(xm) == v0
                xm(boundaryInds)==1
        cvx_end
        
        cvx_begin
            cvx_precision best;
            cvx_solver mosek;
            variable xm(size(sigmavec));
            z4 = xm'*A*xm/2 - b'*xm;
            minimize z4
            subject to
                xm >= 0
                xm <= 1
                sum(xm) == v0
                xm(boundaryInds)==1
        cvx_end
        %}
        
        cvx_begin
            cvx_precision best;
            cvx_solver mosek;
            variable sigmam(size(sigma)) diagonal;
            sigmavecm = diag(sigmam);
            zobj1 = pow_pos(norm(phihat-M*phi,'fro'),2) + trace((D'*sigmam*D*phi-J)*Z') + w*(sum(diag(sigmam))-v0)...
                    + rho1/2*pow_pos(norm(D'*sigmam*D*phi-J,'fro'),2) + rho2/2*pow_abs(sum(diag(sigmam))-v0,2) + ...
                    rho1/2*pow_pos(norm(sigmavecm-diag(sigmaT)),2) + z'*diag(sigmaT-sigmam);
            minimize zobj1
            subject to
                sigmam == diag(diag(sigmam));
                sigmavecm(boundaryEdgeInds)==1
        cvx_end
        zobj2 = pow_pos(norm(phihat-M*phi,'fro'),2) + trace((D'*sigma*D*phi-J)*Z') + w*(sum(diag(sigma))-v0)...
        + rho1/2*pow_pos(norm(D'*sigma*D*phi-J,'fro'),2) + rho2/2*pow_abs(sum(diag(sigma))-v0,2) + ...
        rho1/2*pow_pos(norm(diag(sigma)-diag(sigmaT)),2) + z'*diag(sigmaT-sigma);
        assert(norm((zobj1-zobj2)/zobj2)<1e-2);
    end
    
    % post opt
    times(counter) = toc(tperiter);
    counter = counter + 1;
end

[sortedC, sortedInds] = sort(sigmavec);
keepinds = sortedInds(1:(numel(conductances0)-v0));
figure; hold all; rotate3d on; 
xlim(BB(:,1)'+[1 -1]*1e-1);ylim(BB(:,2)'+[1 -1]*1e-1);zlim(BB(:,3)'+[1 -1]*1e-1);
scatter3(HMesh.edgeCenters(keepinds,1),HMesh.edgeCenters(keepinds,2),HMesh.edgeCenters(keepinds,3),30,ones(numel(keepinds),1),'green','filled');
thresholdedConductances = ones(numel(conductances0),1); thresholdedConductances(keepinds)=0; showinds = find((thresholdedConductances-conductancesGT)~=0);
scatter3(HMesh.edgeCenters(showinds,1),HMesh.edgeCenters(showinds,2),HMesh.edgeCenters(showinds,3),5,thresholdedConductances(showinds)-conductancesGT(showinds));
title('Diff with GT'); xlabel('Green is the resulting reconstruction. Yellow is the missing parts. Blue are the extra parts.');
ptc2 = patch('Faces',HollowHMesh.F2V(HollowHMesh.isBoundaryFace,:),'Vertices',HollowHMesh.V2P,'FaceColor','green','EdgeColor','black'); alpha(ptc2,.05);
%scatter3(HMesh.edgeCenters(knownInds,1),HMesh.edgeCenters(knownInds,2),HMesh.edgeCenters(knownInds,3),5,thresholdedConductances(knownInds)-conductancesGT(knownInds));

















