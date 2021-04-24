function [Op, M, Bnn, clampBasis] = FFOp2D(meshData, Tij, neumann)

if nargin < 3
    neumann = false;
end

nf = meshData.nf;
nv = meshData.nv;
nb = meshData.nb;

%% Construct phase field operators
faceEdgeVecs = reshape(meshData.verts(meshData.faceOrientedEdges(:, 2), :) - ...
                       meshData.verts(meshData.faceOrientedEdges(:, 1), :), nf, 3, 3);
faceEdgeLengths = vecnorm(faceEdgeVecs, 2, 3);
faceEdgeTangents = faceEdgeVecs ./ faceEdgeLengths;
faceEdgeTangents = faceEdgeTangents(:, :, 1) + 1i .* faceEdgeTangents(:, :, 2);
faceEdgeNormals = 1i .* faceEdgeTangents;
faceEdgeNormals = permute(cat(3, real(faceEdgeNormals), imag(faceEdgeNormals)), [3 1 2]);

% Gradient operator
faceEdgeAltitudes = 2 * meshData.areas ./ faceEdgeLengths;
hatGradients = faceEdgeNormals ./ reshape(faceEdgeAltitudes, 1, nf, 3);
gI = repmat((1:2*nf).', 1, 3);
gJ = repelem(meshData.faces(:, [3 1 2]), 2, 1); % Verts opposite edges
G = sparse(gI(:), gJ(:), hatGradients(:), 2 * nf, nv);

% Divergence operator
dI = repmat((1:2:2*nf).' + [0 0 1 1], 1, 1, 3);
dJ = reshape(4 * (meshData.faces(:, [3 1 2]) - 1), nf, 1, 3) + [1 2 3 4];
dij = repmat(permute(hatGradients, [2 1 3]), 1, 2, 1);
D = sparse(dI(:), dJ(:), dij(:), 2*nf, 4*nv);

% Face area weights
A = sparse((1:2*nf).', (1:2*nf).', repelem(meshData.areas, 2, 1), 2*nf, 2*nf);

% Vertex weights
M4 = kron(meshData.star0lump, speye(4));

if neumann
    %% Impose 0-Neumann Boundary Conditions (weakly)
    bdryVtxNormals = meshData.vtxNormals(:, 1:2).';
    bdryVtxTangents = [0 -1; 1 0] * bdryVtxNormals;

    nt = batchop('mult', reshape(bdryVtxNormals, 2, 1, nb), reshape(bdryVtxTangents, 1, 2, nb));
    Bij = [reshape(nt, 1, 4, nb); reshape(multitransp(nt), 1, 4, nb)];
    
    nn = reshape(batchop('mult', reshape(bdryVtxNormals, 2, 1, nb), reshape(bdryVtxNormals, 1, 2, nb)), 4, nb);
    Bnn = sparse(repmat((1:nb), 4, 1), (4 * meshData.bdryIdx).' + (-3:0).', nn, nb, 4 * nv);
    Bnn = Bnn * (D' * A * G);

    BT = batchop('mult', Bij, Tij(:,:,meshData.bdryIdx));
    BTB = batchop('mult', BT, Bij, 'N', 'T');
    BTBi = batchop('pinv', BTB, 2);
    Tij(:,:,meshData.bdryIdx) = Tij(:,:,meshData.bdryIdx) - batchop('mult', batchop('mult', BT, BTBi, 'T', 'N'), BT, 'N', 'N');
end

vertBaseIdx = reshape(4 * (0:nv - 1), 1, 1, nv);
T_I = repmat(vertBaseIdx + (1:4).', 1, 4, 1, 1);
T_J = repmat(vertBaseIdx + (1:4), 4, 1, 1, 1);
T = sparse(T_I(:), T_J(:), Tij(:), 4 * nv, 4 * nv);
TM = T / M4;

if ~neumann
    %% Natural BCs - only take interior nodes
    D = D(:, 4 * meshData.intIdx + (-3:0).');
    TM = TM(4 * meshData.intIdx + (-3:0).', 4 * meshData.intIdx + (-3:0).');
end

%% Form orthotropic thin plate operator
DAG = D' * A * G;
Op = DAG' * TM * DAG;
M = meshData.star0lump;

if nargout > 3
    %% Impose 0-Neumann boundary conditions (clamping)
    faceEdgeLengths = reshape(faceEdgeLengths.', 3, 1, nf);
    B = abs(meshData.d1.');
    B = B(meshData.bdryEdgeIdx, :);
    [bdryFaceIdx, ~] = find(B.');
    bdryFaces = meshData.faces(bdryFaceIdx, :);
    [bdryFaceEdgeIdx, ~] = find((meshData.bdryEdgeIdx == meshData.face2edge(bdryFaceIdx, :)).');
    bdryFaceEdgeIdx = reshape(bdryFaceEdgeIdx, [], 1);
    bdryFaceShift = mod(bdryFaceEdgeIdx + (-2:0), 3) + 1;
    bdryFaceShift = sub2ind(size(bdryFaces), repmat((1:nb).', 1, 3), bdryFaceShift);
    bdryFaces = reshape(bdryFaces(bdryFaceShift), nb, 3);
    bdryFaceVertAngles = meshData.faceVertAngles(bdryFaceIdx, :);
    bdryFaceVertAngles = reshape(bdryFaceVertAngles(bdryFaceShift), nb, 3);
    bdryFaceLengths = squeeze(faceEdgeLengths(:, :, bdryFaceIdx)).';
    bdryFaceLengths = reshape(bdryFaceLengths(bdryFaceShift), nb, 3);

    % vk = bdryFaces(:, 1);
    bdryBasis = sparse(bdryFaces(:, [1 1]), bdryFaces(:, [3 2]), cos(bdryFaceVertAngles(:, [2 3])) .* bdryFaceLengths(:, [1 3]) ./ bdryFaceLengths(:, 2), nv, nv);
    % allButVk = setdiff((1:nv).', vk);
    % bdryBasis = bdryBasis(:, allButVk);
    bdryNeighborCount = sum(bdryBasis ~= 0, 2);
    bdryClampedMask = bdryNeighborCount == 2;
    bdryBasis(~bdryClampedMask, :) = 0;

    %% Critical points at singularities
    ang = angle(z4);
    angDiff = mod(pi + meshData.d1.' * ang, 2 * pi) - pi;
    holonomy = meshData.d0.' * angDiff;
    singIdx = meshData.intIdx(abs(holonomy(meshData.intIdx)) > pi);
    regIdx = setdiff(1:nv, singIdx);
    singBasis = meshData.L;
    singBasis(:, regIdx) = 0;
    singBasis = singBasis ~= 0;
    singNbdIdx = setdiff(find(any(singBasis, 2)), singIdx);
    clampedIdx = unique([find(bdryClampedMask); singNbdIdx]);
    unclampedIdx = setdiff(1:nv, clampedIdx);

    unclampedIdent = speye(nv);
    unclampedIdent(clampedIdx, :) = 0;
    clampBasis = bdryBasis + singBasis + unclampedIdent;
    clampBasis = clampBasis(:, unclampedIdx);
end

end