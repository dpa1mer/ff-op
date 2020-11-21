function [V, lambda] = PhaseField2D(verts, faces, ellipticity, singPenalty)

if nargin < 3 || isempty(ellipticity)
    ellipticity = 0.01;
end

if nargin < 4 || isempty(singPenalty)
    singPenalty = 0;
end

%% Load mesh. Note: triangles must be oriented consistently
faces = flip_ears(verts, faces);
tri = triangulation(faces, verts);
bdryIdx = freeBoundary(tri);
nb = size(bdryIdx, 1);
nv = size(verts, 1);
nf = size(faces, 1);
intIdx = setdiff(1:nv, bdryIdx);

faceCenters = (1/3) * squeeze(sum(reshape(verts(faces, :), nf, 3, 3), 2));

faceOrientedEdges = reshape(faces(:, [1 2; 2 3; 3 1]), [], 2);
faceOrientation = reshape(sign(faceOrientedEdges(:, 2) - faceOrientedEdges(:, 1)), [], 3);
faceEdges = sort(faceOrientedEdges, 2);
[edges, ~, face2edge] = unique(faceEdges, 'rows');
face2edge = reshape(face2edge, [], 3);
ne = size(edges, 1);

%% Construct Hodge Stars
vertOppEdges = face2edge(:, [2 3 1]);
eij = verts(faces(:, [3 1 2]), :) - verts(faces, :);
eik = verts(faces(:, [2 3 1]), :) - verts(faces, :);
vertAngles = reshape(acos((dot(eij, eik, 2) ./ vecnorm(eij, 2, 2)) ./ vecnorm(eik, 2, 2)), nf, 3);
oppositeCotans = cot(vertAngles);

Lij = oppositeCotans / 2;
L = sparse(edges(vertOppEdges, 1), edges(vertOppEdges, 2), Lij, nv, nv);
star1 = spdiags(nonzeros(L'), 0, ne, ne);
star1dual = spdiags(diag(star1).^(-1), 0, ne, ne);

areas = 0.5 * vecnorm(cross(verts(faces(:, 1), :) - verts(faces(:, 2), :), verts(faces(:, 3), :) - verts(faces(:, 2), :)), 2, 2);
star2 = spdiags(areas, 0, nf, nf);

star0 = massmatrix(verts, faces, 'full');%sparse(faces, faces, repmat(areas / 3, 1, 3), nv, nv);
star0lump = massmatrix(verts, faces, 'barycentric');

%% Construct Laplacians
d0 = sparse(repmat((1:ne).', 1, 2), edges, repmat([-1 1], ne, 1), ne, nv);
d1 = sparse(repmat((1:nf).', 1, 3), face2edge, faceOrientation, nf, ne);
L = d0.' * star1 * d0;
Ldual = d1 * star1dual * d1.';

%% Construct boundary
[bdryForward, bdryEdgeIdxF] = ismember(bdryIdx, edges, 'rows');
[bdryBackward, bdryEdgeIdxB] = ismember(fliplr(bdryIdx), edges, 'rows');
bdryEdgeIdx = bdryEdgeIdxF + bdryEdgeIdxB;
Bij = abs(d1.');
Bij = Bij(bdryEdgeIdx, :);
% B = B * star2;

%% Compute cross field boundary conditions
bdryEdgeVecs = verts(bdryIdx(:, 2), :) - verts(bdryIdx(:, 1), :);
bdryNormals = [-bdryEdgeVecs(:, 2), bdryEdgeVecs(:, 1), bdryEdgeVecs(:, 3)];
bdryNormals = bdryNormals ./ vecnorm(bdryNormals, 2, 2);
bdryEdgeCenters = 0.5 * (verts(bdryIdx(:, 2), :) + verts(bdryIdx(:, 1), :));
zb = (bdryNormals(:, 1) + 1i.*bdryNormals(:, 2)).^4;

%% Compute cross field by MBO
z4 = [Ldual Bij.'; Bij zeros(size(Bij, 1))] \ [zeros(nf, 1); zb];
z4 = z4(1:nf);
z4 = z4 ./ abs(z4);

lambda = eigs(Ldual, star2, 2, 'smallestabs');
tau = 0.5/lambda(2);
Ambo = star2 + tau * Ldual;
for j = 1:10
    z4 = [Ambo, Bij.'; Bij zeros(size(Bij, 1))] \ [star2 * z4; zb];
    z4 = z4(1:nf);
    crossNorm = abs(z4);
    z4 = z4 ./ crossNorm;
end

z = (z4.^(1/4)) .* [1 1i -1 -1i];
crossField = cat(3, real(z), imag(z), zeros(size(z)));

% figure; FancyQuiver(repmat(faceCenters, 4, 1), reshape((crossNorm + 1).^(-1) .* crossField, nf * 4, 3), inferno, 0); view(2); axis image off;

%% Construct phase field operators
faceEdgeVecs = reshape(verts(faceOrientedEdges(:, 2), :) - verts(faceOrientedEdges(:, 1), :), nf, 3, 3);
faceEdgeLengths = vecnorm(faceEdgeVecs, 2, 3);
faceEdgeTangents = faceEdgeVecs ./ faceEdgeLengths;
faceEdgeTangents = faceEdgeTangents(:, :, 1) + 1i .* faceEdgeTangents(:, :, 2);
faceEdgeNormals = 1i .* faceEdgeTangents;
% faceEdgeTangents = permute(cat(3, real(faceEdgeTangents), imag(faceEdgeTangents)), [3 1 2]);
faceEdgeNormals = permute(cat(3, real(faceEdgeNormals), imag(faceEdgeNormals)), [3 1 2]);

% Gradient operator
faceEdgeAltitudes = 2 * areas ./ faceEdgeLengths;
hatGradients = faceEdgeNormals ./ reshape(faceEdgeAltitudes, 1, nf, 3);
gI = repmat((1:2*nf).', 1, 3);
gJ = repelem(faces(:, [3 1 2]), 2, 1); % Verts opposite edges
G = sparse(gI(:), gJ(:), hatGradients(:), 2 * nf, nv);

% Divergence operator
dI = repmat((1:2:2*nf).' + [0 0 1 1], 1, 1, 3);
dJ = reshape(4 * (faces(:, [3 1 2]) - 1), nf, 1, 3) + [1 2 3 4];
dij = repmat(permute(hatGradients, [2 1 3]), 1, 2, 1);
D = sparse(dI(:), dJ(:), dij(:), 2*nf, 4*nv);

% Face area weights
A = sparse((1:2*nf).', (1:2*nf).', repelem(areas, 2, 1), 2*nf, 2*nf);

% Vertex weights
M4 = kron(star0lump, speye(4));

% Cross field tensor
u = [real(z(:, 1)), imag(z(:, 1))].';
v = [real(z(:, 2)), imag(z(:, 2))].';
u2 = reshape(reshape(u, 2, 1, nf) .* reshape(u, 1, 2, nf), 4, 1, nf);
u4 = u2 .* reshape(u2, 1, 4, nf);
v2 = reshape(reshape(v, 2, 1, nf) .* reshape(v, 1, 2, nf), 4, 1, nf);
v4 = v2 .* reshape(v2, 1, 4, nf);
Tij = eye(4) - (1 - ellipticity) * (u4 + v4);
% Tij = reshape((crossNorm).^(-1), 1, 1, nf) .* Tij;
areaWeightedT = repmat(reshape(areas, 1, 1, nf) .* Tij, 1, 1, 1, 3) / 3;
faceBaseIdx = 4 * (reshape(faces, 1, 1, nf, 3) - 1);
MT_I = repmat(faceBaseIdx + [1;2;3;4], 1, 4, 1, 1);
MT_J = repmat(faceBaseIdx + [1 2 3 4], 4, 1, 1, 1);
MT = sparse(MT_I(:), MT_J(:), areaWeightedT(:), 4 * nv, 4 * nv);
T = M4 \ MT;


%% Impose 0-Neumann Boundary Conditions (weakly)
bdryVtxIdx = unique(bdryIdx);
bdryVtxNormals = accumarray([repmat(bdryIdx(:, 1), 3, 1) repelem((1:3).', nb)], bdryNormals(:), [nv 3]) ...
               + accumarray([repmat(bdryIdx(:, 2), 3, 1) repelem((1:3).', nb)], bdryNormals(:), [nv 3]);
bdryVtxNormals = bdryVtxNormals(bdryVtxIdx, 1:2).';
bdryVtxTangents = [0 -1; 1 0] * bdryVtxNormals;

nt = batchop('mult', reshape(bdryVtxNormals, 2, 1, nb), reshape(bdryVtxTangents, 1, 2, nb));
Bij = [reshape(nt, 1, 4, nb); reshape(multitransp(nt), 1, 4, nb)];

Bi = repmat(reshape(1:2*nb, 2, 1, nb), 1, 4, 1);
Bj = repmat(reshape((4 * bdryVtxIdx + (-3:0)).', 1, 4, nb), 2, 1, 1);
B = sparse(Bi(:), Bj(:), Bij(:), 2 * nb, 4 * nv);

Tadjusted = T - T * (B' * ((B * T * B') \ (B * T)));
TM = Tadjusted / M4;

%% Penalize gradient at singularities
S = G' * (repelem((crossNorm + 0.1).^(-1), 2, 1) .* A) * G;

%% Compute eigenfunctions
DAG = D' * A * G;
O = DAG' * TM * DAG + singPenalty * S;
M = star0lump;
[V, lambda] = eigs(O + 1e-6 * M, M, 200, 'smallestabs');%, 'IsSymmetricDefinite', true);
figure;
% [cmin, cmax] = bounds(V(:, 1:64), 'all');
for k = 1:64
    subplot(8, 8, k);
    trisurf(faces, verts(:, 1), verts(:, 2), verts(:, 3), V(:, k), 'EdgeColor', 'none'); view(2); axis image off; shading interp; colormap viridis; %caxis([cmin cmax]);
end



end