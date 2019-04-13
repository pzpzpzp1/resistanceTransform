
function HMesh = loadHexMesh(V2P, H2V)

HMesh.H2V = H2V;
HMesh.V2P = V2P;
HMesh.nverts = size(V2P,1);
HMesh.nhexes = size(H2V,1);
ncells = size(H2V,1);
npoints = size(V2P,1);
cells = H2V;



%% use vtk to see singular verts.
%npoints x ncells
indicator = sparse(repmat([1:ncells]',8,1),cells(:),ones(8*ncells,1))';
vertdegs = sum(indicator')';

HMesh.hdegv = vertdegs;

%{
citer = randi(ncells);
pts = points(cells(citer,:)',:);
scatter3(pts(:,1), pts(:,2), pts(:,3));
edgeinds = [1 2 1 4 1 5 2 3 3 4 5 6 6 7 7 8 5 8 2 6 3 7 4 8];
hold all;
for i = 1:(numel(edgeinds)/2)
    v1 = edgeinds((i-1)*2+1);
    v2 = edgeinds((i-1)*2+2);
    v = [pts(v1,:); pts(v2,:)];
    plot3(v(:,1),v(:,2),v(:,3));
    %pause;
end
pause;
%}

edgeinds = [1 2 1 4 1 5 2 3 3 4 5 6 6 7 7 8 5 8 2 6 3 7 4 8];
edgest = reshape(cells(:,edgeinds)', 2, 12*ncells)';
edges = unique(sort(edgest,2),'rows');
sharedHexesPerEdge = indicator(edges(:,1),:) .* indicator(edges(:,2),:);
edgedegs = sum(sharedHexesPerEdge')';
edgeCenters = (HMesh.V2P(edges(:,1),:)+HMesh.V2P(edges(:,2),:))/2;
edgeLengths = norms(HMesh.V2P(edges(:,1),:)-HMesh.V2P(edges(:,2),:),2,2);

HMesh.edgeCenters = edgeCenters;
HMesh.edgeLengths = edgeLengths;
HMesh.hdege = edgedegs;
HMesh.E2V = edges;
HMesh.VV2E = sparse(edges(:,1),edges(:,2),1:size(edges,1),npoints,npoints); HMesh.VV2E = HMesh.VV2E + HMesh.VV2E';
HMesh.nedges = size(edges,1);

faceinds = [1 2 3 4 1 2 6 5 1 4 8 5 7 6 5 8 7 6 2 3 7 8 4 3];
facest = reshape(cells(:,faceinds)', 4, [])';
facesSorted = sort(facest,2); 
[C, I1, I2] = unique(facesSorted,'rows');
faces = facest(I1,:);
sharedHexesPerFace = indicator(faces(:,1),:) .* indicator(faces(:,2),:) .* indicator(faces(:,3),:) .* indicator(faces(:,4),:);
facedegs = sum(sharedHexesPerFace')';
isBoundaryFace = (facedegs==1);  

HMesh.p_edgeinds = edgeinds;
HMesh.p_faceinds = faceinds;
HMesh.F2HM = sharedHexesPerFace;

HMesh.F2V = faces;
HMesh.nfaces = size(faces,1);
HMesh.FaceCenters = (HMesh.V2P(HMesh.F2V(:,1),:) + HMesh.V2P(HMesh.F2V(:,2),:) + HMesh.V2P(HMesh.F2V(:,3),:) + HMesh.V2P(HMesh.F2V(:,4),:))/4;
HMesh.isBoundaryFace = isBoundaryFace;

bfaces = faces(isBoundaryFace,:);
boundaryverts = unique(bfaces(:));
isBoundaryVerts = sparse(boundaryverts,ones(numel(boundaryverts),1), ones(numel(boundaryverts),1),npoints,1)==1;

HMesh.isBoundaryVerts = isBoundaryVerts;

bedges = sort(reshape(bfaces(:,[1 2 2 3 3 4 4 1])', 2, [])',2);
bedges = unique(bedges,'rows');
nbedges = size(bedges,1);

sprsE = sparse(repmat(1:size(edges,1),1,2),edges(:),ones(numel(edges),1),HMesh.nedges,npoints);
sprsBE = sparse(repmat([1:size(bedges,1)],1,2),bedges(:),ones(numel(bedges),1),nbedges,npoints);

% pairwise comparision edges x boundaryedges to find indices.
BEInd = (sprsE*sprsBE')==2;
isBoundaryEdge = sum(BEInd')';

HMesh.isBoundaryEdge = isBoundaryEdge==1;

% construct incidence matrix
rectEdges = HMesh.E2V;
nE = HMesh.nedges;
nV = HMesh.nverts;
nC = HMesh.nhexes;
rectIncidenceMatrix = sparse([1:nE 1:nE]', [rectEdges(:,1);rectEdges(:,2)], [ones(nE,1) -ones(nE,1)], nE, nV);
gradientOp = sparse([1:nE 1:nE]', [rectEdges(:,1);rectEdges(:,2)], [ones(nE,1)./HMesh.edgeLengths -ones(nE,1)./HMesh.edgeLengths], nE, nV);
HMesh.incidenceMatrix = rectIncidenceMatrix;
HMesh.gradientOp = gradientOp;



if false
    isSingularVerts = (isBoundaryVerts & vertdegs ~= 4) | (~isBoundaryVerts & vertdegs ~= 8);
    isSingularEdge = (isBoundaryEdge & edgedegs ~= 2) | (~isBoundaryEdge & edgedegs ~= 4);
    isInteriorDeg3SingularEdge = (~isBoundaryEdge & edgedegs == 3);
    isInteriorDeg5SingularEdge = (~isBoundaryEdge & edgedegs == 5);
    isSingularBoundaryVert = isBoundaryVerts & vertdegs ~= 4;

    HMesh.isSingularEdge = isSingularEdge;
    HMesh.isInteriorDeg3SingularEdge = isInteriorDeg3SingularEdge;
    HMesh.isInteriorDeg5SingularEdge = isInteriorDeg5SingularEdge;


    %% try to build graph of singularity structure.
    % singular vert to vert
    sv2v = edges(find(isSingularEdge), :);
    sVertsAll = unique(sv2v(:));
    sv2vM = sparse(sv2v(:,1), sv2v(:,2), ones(size(sv2v,1),1), npoints, npoints);
    % guaranteed all singular meta vertices.
    sv2vMSym = ((sv2vM + sv2vM') ~= 0);
    sverts = find(sum(sv2vM + sv2vM')~=2 & sum(sv2vM + sv2vM')~=0);

    HMesh.nsingularverts = numel(sv2v);
    HMesh.nsingularMetaverts = numel(sverts);
    HMesh.isSingularVertex = sparse(ones(numel(sVertsAll),1),sVertsAll,ones(numel(sVertsAll),1),1,HMesh.nverts);
    HMesh.isSingularMetaVertex = sparse(ones(numel(sverts),1),sverts,ones(numel(sverts),1),1,HMesh.nverts);
end
%% find corners. useful if you want the singular vertex types. but this segment is copy pasted and may not work.
if false
MetaVertices = cell(numel(sverts),1);
for i = 1:numel(sverts)
    mvert.vind = sverts(i);
    mvert.isBoundary = HMesh.isBoundaryVerts(sverts(i));
    adjHexes = find(sum(HMesh.H2V == sverts(i),2));
    oneRingVerts = find(HMesh.VV2E(sverts(i),:));
    
    corners = {}; pos = 1;
    for j = 1:numel(adjHexes)
        hexnum = adjHexes(j);
        %hexE2V = reshape(HMesh.H2V(hexnum, edgeinds),2,[])';
        hexEdges = find(sharedHexesPerEdge(:,hexnum));
        hexE2V = HMesh.E2V(hexEdges,:);
        hexIsSingularEdge = full(HMesh.isSingularEdge(hexEdges));
        hexEdgeDeg = full(HMesh.hdege(hexEdges));
        hexEdgeIsBoundary = full(HMesh.isBoundaryEdge(hexEdges));
        threeEdgeInds = find(sum(hexE2V==sverts(i),2));
            
        threeEdges2V = hexE2V(threeEdgeInds,:);
        if(all(hexIsSingularEdge(threeEdgeInds)))
            % check right handedness by geometry.
            v1 = HMesh.V2P(threeEdges2V(1,2),:) - HMesh.V2P(threeEdges2V(1,1),:);
            v2 = HMesh.V2P(threeEdges2V(2,2),:) - HMesh.V2P(threeEdges2V(1,1),:);
            v3 = HMesh.V2P(threeEdges2V(3,2),:) - HMesh.V2P(threeEdges2V(1,1),:);
            if(dot(cross(v1,v2),v3)<0)
                threeEdgeInds = threeEdgeInds([2 1 3],:);
            end
            
            threeEdges2V = hexE2V(threeEdgeInds,:);
            %threeEdges = hexEdges(threeEdgeInds);
            %threeEdgesIsSingular = hexIsSingularEdge(threeEdgeInds);
            threeEdgesDeg = hexEdgeDeg(threeEdgeInds);
            threeEdgesIsBoundary = hexEdgeIsBoundary(threeEdgeInds);

            corners{pos}.threeEdges2V = threeEdges2V;
            corners{pos}.threeEdgesDeg = threeEdgesDeg;
            corners{pos}.threeEdgesIsBoundary = threeEdgesIsBoundary;
            pos = pos + 1;
        end
    end
    mvert.corners = corners;
    MetaVertices{i}=mvert;
end
HMesh.MetaVertices = MetaVertices;
end

%{
% find meta-adjacency matrix
sv2vMSymT = sv2vMSym; % this one will be butchered in the process :/
MV2MV = sparse(zeros(1)); MV2MV(numel(sverts),numel(sverts))=0;
for i = 1:numel(sverts)
    [num2str(i) '/' num2str(numel(sverts))]
    directions = find(sv2vMSymT(sverts(i),:));
    startloc = i;
    for j = 1:numel(directions)
        prevloc = sverts(i);
        dir = directions(j);
        %sv2vMSymT(sverts(i),dir)=0; sv2vMSymT(dir,sverts(i))=0;
        % keep pursuing this dir until we reach another meta-vertex.
        propdirs = find(sv2vMSymT(dir,:));
        while(numel(propdirs)==2)
            nextloc = ARemoveB(propdirs,prevloc);
            %sv2vMSymT(dir, nextloc)=0; sv2vMSymT(nextloc, dir)=0;
            prevloc = dir;
            dir = nextloc;
            propdirs = find(sv2vMSymT(dir,:));
        end
        endloc = find(dir==sverts);
        assert(numel(endloc)==1);
        edgeind = find(sum([prevloc dir] == edges,2)==2 | sum([dir prevloc] == edges,2)==2);
        edeg = edgedegs(edgeind);
        
        MV2MV(startloc,endloc)=edeg;
        MV2MV(endloc,startloc)=edeg;
    end
end
[mv1,mv2]=find(MV2MV);
SV2SV = sparse(sverts(mv1), sverts(mv2), MV2MV(find(MV2MV)));
[sv1, sv2] = find(SV2SV);

%% load triangulations representing singular vertices
DEGREEMAX  = 10;
Z = LoadTriangulations('D2_BE6', 'S2', DEGREEMAX);

%% map singular vertices to triangle meshes
SV2T = cell(size(sverts));
for i = 1:numel(sverts)
    sv = sverts(i);
    degrees = nonzeros(SV2SV(sv,:))'; %SV2SV is symmetric already.
    %adjSingularVerts = find(SV2SV(sv,:));
    isBVert = full(isBoundaryVerts(sv));
    
    adjEdges = nonzeros(Adjacency(sv,:)')';
    totalDegrees = full(edgedegs(adjEdges) + isBoundaryEdge(adjEdges))'; % hdeg to fdeg
    
    signature = DegreesToSignature(totalDegrees, DEGREEMAX);
    if(isBVert)
        match = find(sum(signature == Z.D2Signatures,2)==size(signature,2));
        assert(numel(match)~=0);
        SV2T{i} = Z.D2{match(1)};
    else
        match = find(sum(signature == Z.S2Signatures,2)==size(signature,2));
        assert(numel(match)~=0);
        SV2T{i} = Z.S2{match(1)};
    end
    
    if(numel(match)~=1)
        error('sucks paul. cant be lazy. That boundary case actually shows up apparently.');
    end
end
%}

if(false)
    
    %Visualize simplified graph.
    %{
    figure; hold on; axis equal;
    for i=1:numel(sv1)
        s = SV2SV(sv1(i), sv2(i));
        assert(s == 3 | s == 5);
        if(s==3); color = 'g'; else; color = 'r'; end;
        vs = points([sv1(i);sv2(i)],:);
        plot3(vs(:,1),vs(:,2),vs(:,3),color);
    end
    scatter3(points(sverts,1),points(sverts,2),points(sverts,3),50,'k.');
    %}
    
    %{
    hold all;
    for i = 1:size(bedges,1)
        if(mod(i,1000)==0)
            pause
        end
        v = points(bedges(i,:)',:);
        plot3(v(:,1), v(:,2), v(:,3),'r','LineWidth',1);
    end
    %}

    singedginds = find(isSingularEdge);
    Ising3 = find(isInteriorDeg3SingularEdge);
    Ising5 = find(isInteriorDeg5SingularEdge);
    beinds = find(isBoundaryEdge);
    hold on; axis equal;
    for i=1:size(singedginds,1)
        v = points(edges(singedginds(i),:)',:);
        plot3(v(:,1), v(:,2), v(:,3),'b','LineWidth',1);
        %pause(.01)
    end
    for i=1:size(Ising3,1)
        v = points(edges(Ising3(i),:)',:);
        plot3(v(:,1), v(:,2), v(:,3),'g','LineWidth',3);
    end
    for i=1:size(Ising5,1)
        v = points(edges(Ising5(i),:)',:);
        plot3(v(:,1), v(:,2), v(:,3),'r','LineWidth',3);
    end
    %scatter3(points(boundaryverts,1),points(boundaryverts,2),points(boundaryverts,3),1,'k.');
    scatter3(points(find(isSingularBoundaryVert),1),points(find(isSingularBoundaryVert),2),points(find(isSingularBoundaryVert),3),100,'b.');
    %scatter3(points(verts,1),points(verts,2),points(verts,3),1,'k.');

    %{
    for i=1:size(beinds,1)
        %if(mod(i,1000)==0)
        %    pause
        %end
        v = points(edges(beinds(i),:)',:);
        plot3(v(:,1), v(:,2), v(:,3),'k','LineWidth',.0001);
    end
    %}
    %{
    hold on;
    for i=1:size(edges,1)
    %     if(mod(i,1000)==0)
    %        pause
    %     end
        v = points(edges(i,:)',:);
        plot3(v(:,1), v(:,2), v(:,3),'k','LineWidth',.0001);
        pause(.0001);
    end
    %}

    %scatter3(points(isSingularVerts,1),points(isSingularVerts,2),points(isSingularVerts,3),'k.');
end

end




