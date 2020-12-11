function [Op, star0lump] = PhaseField3D(meshData, ellipticity)

if nargin < 2
    ellipticity = 0.01;
end

%% Load mesh.
assert(~meshData.isGrid);
L = meshData.L;
star0 = meshData.M;

verts = meshData.verts;
tets = meshData.tets;
nt = meshData.nt;
nv = meshData.nv;
bdryIdx = meshData.bdryIdx;
nb = length(bdryIdx);

%% Compute frame field by MBO

q = MBO(meshData, RayMBO);
octa = OctaMBO;
q = octa.proj(q);
% frames = Coeff2Frames(q);

% Convert frames to tensors
frameTensors = reshape(Monomial2Tensor(Sph024ToMonomial(Octa2Odeco(q))), 9, 9, nv);
Tij = eye(9) - (1 - ellipticity) * (frameTensors ./ 24);

% VisualizeFrameField(meshData, frames);

%% Construct phase field operators
volumes = TetVolumes(verts, tets);

tetOppFaces = reshape(tets(:, flipud(nchoosek(1:4, 3))), nt, 4, 3);
oppFaceNormals = cross(verts(tetOppFaces(:, :, 2), :) - verts(tetOppFaces(:, :, 1), :), ...
                       verts(tetOppFaces(:, :, 3), :) - verts(tetOppFaces(:, :, 1), :), 2);
oppFaceCenters = squeeze(sum(reshape(verts(tetOppFaces, :), nt, 4, 3, 3), 3) ./ 3);
oppFaceAreas = vecnorm(oppFaceNormals, 2, 2) ./ 2;
oppFaceNormals = oppFaceNormals ./ oppFaceAreas;
normalSigns = sign(dot(verts(tets, :) - verts(tetOppFaces(:, :, 1), :), oppFaceNormals, 2));
oppFaceNormals = normalSigns .* oppFaceNormals;
oppFaceNormals = reshape(oppFaceNormals, nt, 4, 3);
oppFaceAltitudes = 3 * reshape(volumes, nt, 1) ./ reshape(oppFaceAreas, nt, 4);

% Gradient operator
hatGradients = permute(oppFaceNormals ./ oppFaceAltitudes, [3 2 1]);
gI = repmat(reshape(1:3*nt, 3, 1, nt), 1, 4, 1);
gJ = repmat(reshape(tets.', 1, 4, nt), 3, 1, 1);
G = sparse(gI(:), gJ(:), hatGradients(:), 3 * nt, nv);

% Divergence operator
dI = repmat((1:3:3*nt).' + [0 0 0 1 1 1 2 2 2], 1, 1, 4);
dJ = reshape(9 * (tets - 1), nt, 1, 4) + (1:9);
dij = repmat(permute(hatGradients, [3 1 2]), 1, 3, 1);
D = sparse(dI(:), dJ(:), dij(:), 3*nt, 9*nv);

% Tet volume weights
V = spdiags(repelem(volumes, 3, 1), 0, 3 * nt, 3 * nt);

% Vertex weights
star0lump = sparse(tets, tets, repmat(volumes ./ 4, 1, 4), nv, nv);
M9 = kron(star0lump, speye(9));

% weightedT = reshape(full(diag(star0lump)), 1, 1, nv) .* Tij;

%% Impose 0-Neumann Boundary Conditions (weakly)
% bdryNormals = reshape(meshData.bdryNormals, 3, 1, nb);
% nn = batchop('mult', bdryNormals, bdryNormals, 'N', 'T');
% nnn = batchop('mult', bdryNormals, reshape(nn, 1, 9, nb));
% multN = [multitransp(bdryNormals) zeros(1, 6, nb);
%          zeros(1, 3, nb) multitransp(bdryNormals) zeros(1, 3, nb);
%          zeros(1, 6, nb) multitransp(bdryNormals)];
% multNt = multitransp(reshape(multN, 9, 3, nb));
% Bij = [multN - nnn; multNt - nnn];
% BT = batchop('mult', Bij, Tij(:,:,bdryIdx));
% BTB = batchop('mult', BT, Bij, 'N', 'T');
% BTBi = batchop('pinv', BTB, 5);
% Tij(:,:,bdryIdx) = Tij(:,:,bdryIdx) - batchop('mult', batchop('mult', BT, BTBi, 'T', 'N'), BT, 'N', 'N');

% Bi = repmat(reshape(1:6*nb, 6, 1, nb), 1, 9, 1);
% Bj = repmat(reshape((9 * bdryIdx + (-8:0)).', 1, 9, nb), 6, 1, 1);
% B = sparse(Bi(:), Bj(:), Bij(:), 6 * nb, 9 * nv);

%% Form frame field principal symbol matrix
vertBaseIdx = reshape(9 * (0:nv - 1), 1, 1, nv);
T_I = repmat(vertBaseIdx + (1:9).', 1, 9, 1, 1);
T_J = repmat(vertBaseIdx + (1:9), 9, 1, 1, 1);
T = sparse(T_I(:), T_J(:), Tij(:), 9 * nv, 9 * nv);
TM = T / M9;

% % Only take interior nodes
% intBlkIdx = (9 * intIdx + (-8:0)).';
% D = D(:, intBlkIdx);
% M9 = M9(intBlkIdx, intBlkIdx);
% MT = MT(intBlkIdx, intBlkIdx);

% %% Impose 0-Neumann boundary conditions
% 
% 
% faceEdgeLengths = reshape(faceEdgeLengths.', 3, 1, nf);
% [bdryFaceIdx, ~] = find(B.');
% bdryFaces = faces(bdryFaceIdx, :);
% [bdryFaceEdgeIdx, ~] = find((bdryEdgeIdx == face2edge(bdryFaceIdx, :)).');
% bdryFaceEdgeIdx = reshape(bdryFaceEdgeIdx, [], 1);
% bdryFaceShift = mod(bdryFaceEdgeIdx + (-2:0), 3) + 1;
% bdryFaceShift = sub2ind(size(bdryFaces), repmat((1:nb).', 1, 3), bdryFaceShift);
% bdryFaces = reshape(bdryFaces(bdryFaceShift), nb, 3);
% bdryFaceVertAngles = vertAngles(bdryFaceIdx, :);
% bdryFaceVertAngles = reshape(bdryFaceVertAngles(bdryFaceShift), nb, 3);
% bdryFaceLengths = squeeze(faceEdgeLengths(:, :, bdryFaceIdx)).';
% bdryFaceLengths = reshape(bdryFaceLengths(bdryFaceShift), nb, 3);
% 
% % vk = bdryFaces(:, 1);
% bdryBasis = sparse(bdryFaces(:, [1 1]), bdryFaces(:, [3 2]), cos(bdryFaceVertAngles(:, [2 3])) .* bdryFaceLengths(:, [1 3]) ./ bdryFaceLengths(:, 2), nv, nv);
% % allButVk = setdiff((1:nv).', vk);
% % bdryBasis = bdryBasis(:, allButVk);
% bdryNeighborCount = sum(bdryBasis ~= 0, 2);
% constrainedNodes = bdryNeighborCount == 2;
% bdryBasis(~constrainedNodes, :) = 0;
% otherIdent = speye(nv);
% otherIdent(constrainedNodes, :) = 0;
% bdryBasis = bdryBasis + otherIdent;
% bdryBasis = bdryBasis(:, ~constrainedNodes);

%% Compute eigenfunctions
DVG = D' * V * G;
Op = DVG' * TM * DVG;
% O = bdryBasis.' * O * bdryBasis;
% M = bdryBasis.' * star0 * bdryBasis;
% [V, lambda] = eigs(Op + star0lump, star0lump, 200, 'smallestabs');%, 'IsSymmetricDefinite', true);
% % V = bdryBasis * V;
% % [W, ~] = eigs(L, star0, 100, 'smallestabs');
% figure;
% % k = 200;
% % [cmin, cmax] = bounds(V(:, 1:64), 'all');
% for k = 1:64
%     subplot(8, 8, k);
%     trisurf(meshData.bdry, V(:, k)); view(3); axis image vis3d off; shading interp; colormap viridis;
% end

end