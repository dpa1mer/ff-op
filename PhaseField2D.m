function PhaseField2D(verts, faces)

%% Load mesh. Note: triangles must be oriented consistently
% [verts, faces] = load_mesh('~/Downloads/clover.obj');
tri = triangulation(faces, verts);
bdryIdx = freeBoundary(tri);
nb = size(bdryIdx, 1);
nv = size(verts, 1);
nf = size(faces, 1);

faceCenters = (1/3) * squeeze(sum(reshape(verts(faces, :), nf, 3, 3), 2));

faceOrientedEdges = reshape(faces(:, [1 2; 2 3; 3 1]), [], 2);
faceOrientation = reshape(sign(faceOrientedEdges(:, 2) - faceOrientedEdges(:, 1)), [], 3);
faceEdges = sort(faceOrientedEdges, 2);
[edges, ~, face2edge] = unique(faceEdges, 'rows');
face2edge = reshape(face2edge, [], 3);
ne = size(edges, 1);

%% Construct Hodge Stars
vertOppEdges = reshape(faces(:, [2 3; 1 3; 1 2]), nf, 3, 2);
eij = verts(faces(:, [3 1 2]), :) - verts(faces, :);
eik = verts(faces(:, [2 3 1]), :) - verts(faces, :);
vertAngles = reshape(acos((dot(eij, eik, 2) ./ vecnorm(eij, 2, 2)) ./ vecnorm(eik, 2, 2)), nf, 3);
oppositeCotans = cot(vertAngles);

Lij = oppositeCotans / 2;
L = sparse(vertOppEdges(:, :, 1), vertOppEdges(:, :, 2), Lij, nv, nv);
star1 = spdiags(nonzeros(L'), 0, ne, ne);
star1dual = spdiags(diag(star1).^(-1), 0, ne, ne);

areas = 0.5 * vecnorm(cross(verts(faces(:, 1), :) - verts(faces(:, 2), :), verts(faces(:, 3), :) - verts(faces(:, 2), :)), 2, 2);
star2 = spdiags(areas, 0, nf, nf);

star0 = sparse(faces, faces, repmat(areas / 3, 1, 3), nv, nv);

%% Construct Laplacians
d0 = sparse(repmat((1:ne).', 1, 2), edges, repmat([-1 1], ne, 1), ne, nv);
d1 = sparse(repmat((1:nf).', 1, 3), face2edge, faceOrientation, nf, ne);
L = d0.' * star1 * d0;
Ldual = d1 * star1dual * d1.';

%% Construct boundary
[bdryForward, bdryEdgeIdxF] = ismember(bdryIdx, edges, 'rows');
[bdryBackward, bdryEdgeIdxB] = ismember(fliplr(bdryIdx), edges, 'rows');
bdryEdgeIdx = bdryEdgeIdxF + bdryEdgeIdxB;
B = abs(d1.');
B = B(bdryEdgeIdx, :);
% B = B * star2;

%% Compute cross field boundary conditions
bdryEdgeVecs = verts(bdryIdx(:, 2), :) - verts(bdryIdx(:, 1), :);
bdryNormals = [-bdryEdgeVecs(:, 2), bdryEdgeVecs(:, 1), bdryEdgeVecs(:, 3)];
bdryNormals = bdryNormals ./ vecnorm(bdryNormals, 2, 2);
bdryEdgeCenters = 0.5 * (verts(bdryIdx(:, 2), :) + verts(bdryIdx(:, 1), :));
zb = (bdryNormals(:, 1) + 1i.*bdryNormals(:, 2)).^4;

%% Compute cross field
z4 = [Ldual B.'; B zeros(size(B, 1))] \ [zeros(nf, 1); zb];
z4 = z4(1:nf);
z4 = z4 ./ abs(z4);

z = (z4.^(1/4)) .* [1 1i -1 -1i];
crossField = cat(3, real(z), imag(z), zeros(size(z)));

figure; FancyQuiver(repmat(faceCenters, 4, 1), reshape(crossField, nf * 4, 3), inferno, 0); view(2); axis image off;

%% Construct phase field operators
faceEdgeVecs = reshape(verts(faceOrientedEdges(:, 2), :) - verts(faceOrientedEdges(:, 1), :), nf, 3, 3);
faceEdgeLengths = vecnorm(faceEdgeVecs, 2, 3);
faceEdgeTangents = faceEdgeVecs ./ faceEdgeLengths;
faceEdgeTangents = faceEdgeTangents(:, :, 1) + 1i .* faceEdgeTangents(:, :, 2);
faceEdgeNormals = 1i .* faceEdgeTangents;

% Express in local frame
faceEdgeTangents = faceEdgeTangents ./ z(:, 1);
faceEdgeNormals = faceEdgeNormals ./ z(:, 1);

faceEdgeLengths = reshape(faceEdgeLengths.', 3, 1, nf);
faceEdgeTangents = permute(cat(3, real(faceEdgeTangents), imag(faceEdgeTangents)), [2 3 1]);
faceEdgeNormals = permute(cat(3, real(faceEdgeNormals), imag(faceEdgeNormals)), [2 3 1]);

% Crouzeix-Raviart Frame-Dirichlet matrix
% Isotropic Part
crOrientation = repmat(reshape(faceOrientation.', 3, 1, nf), 2, 1, 1);
crTGrad = reshape(reshape(faceEdgeTangents, 3, 2, 1, nf) .* reshape(faceEdgeNormals, 3, 1, 2, nf), 3, 4, nf);
crNGrad = reshape(reshape(faceEdgeNormals, 3, 2, 1, nf) .* reshape(faceEdgeNormals, 3, 1, 2, nf), 3, 4, nf);
crTNGrad = [crTGrad; crNGrad];
crTNGrad = crTNGrad .* crOrientation;

% Anisotropic Part
crTGradContracted = faceEdgeTangents .* faceEdgeNormals;
crNGradContracted = faceEdgeNormals .* faceEdgeNormals;
crTNGradContracted = [crTGradContracted; crNGradContracted];
crTNGradContracted = crTNGradContracted .* crOrientation;

areas = reshape(areas, 1, 1, nf);
crTNxTN = (batchop('mult', crTNGrad, crTNGrad, 'N', 'T') ...
         - batchop('mult', crTNGradContracted, crTNGradContracted, 'N', 'T')) ./ areas;

faceEdgeTNIdx = [2 * face2edge - 1, 2 * face2edge].';
crI = repmat(reshape(faceEdgeTNIdx, 6, 1, nf), 1, 6, 1);
crJ = repmat(reshape(faceEdgeTNIdx, 1, 6, nf), 6, 1, 1);
CR = sparse(crI(:), crJ(:), crTNxTN(:), 2 * ne, 2 * ne);

% Crouzeix-Raviart mass matrix
crMij = (areas ./ 3) .* repmat((faceEdgeLengths).^(-2), 2, 1, 1) .* eye(6);
Mcr = sparse(crI(:), crJ(:), crMij(:), 2 * ne, 2 * ne);

% Differential matrix
hatGrad = (faceEdgeLengths ./ 6) .* faceEdgeNormals;
hatGrad = hatGrad([2 3 1], :, :);

crTNForm = [faceEdgeTangents ./ faceEdgeLengths; faceEdgeNormals ./ faceEdgeLengths] .* crOrientation;

crDij = batchop('mult', crTNForm, hatGrad, 'N', 'T');
crDI = repmat(reshape(faceEdgeTNIdx, 6, 1, nf), 1, 3, 1);
crDJ = repmat(reshape(faces.', 1, 3, nf), 6, 1, 1);
Dcr = sparse(crDI(:), crDJ(:), crDij(:), 2 * ne, nv);

%% Impose 0-Neumann boundary conditions

[bdryFaceIdx, ~] = find(B.');
bdryFaces = faces(bdryFaceIdx, :);
[bdryFaceEdgeIdx, ~] = find((bdryEdgeIdx == face2edge(bdryFaceIdx, :)).');
bdryFaceEdgeIdx = reshape(bdryFaceEdgeIdx, [], 1);
bdryFaceShift = mod(bdryFaceEdgeIdx + (-2:0), 3) + 1;
bdryFaceShift = sub2ind(size(bdryFaces), repmat((1:nb).', 1, 3), bdryFaceShift);
bdryFaces = reshape(bdryFaces(bdryFaceShift), nb, 3);
bdryFaceVertAngles = vertAngles(bdryFaceIdx, :);
bdryFaceVertAngles = reshape(bdryFaceVertAngles(bdryFaceShift), nb, 3);
bdryFaceLengths = squeeze(faceEdgeLengths(:, :, bdryFaceIdx)).';
bdryFaceLengths = reshape(bdryFaceLengths(bdryFaceShift), nb, 3);

bdryBasis = speye(nv);
vk = bdryFaces(:, 1);
bdryBasis = bdryBasis + sparse(bdryFaces(:, [1 1]), bdryFaces(:, 2:3), cos(bdryFaceVertAngles(:, 2:3)) .* bdryFaceLengths(:, [1 3]) ./ bdryFaceLengths(:, 2), nv, nv);
allButVk = setdiff((1:nv).', vk);
bdryBasis = bdryBasis(:, allButVk);

%% Compute eigenfunctions
A = (Dcr.' * (Mcr \ (CR * (Mcr \ Dcr))));
A = bdryBasis.' * A * bdryBasis;
M = bdryBasis.' * star0 * bdryBasis;
[V, lambda] = eigs(A, M, 100, 'smallestabs');%, 'IsSymmetricDefinite', true);
V = bdryBasis * V;
% [W, ~] = eigs(L, star0, 100, 'smallestabs');
figure;
% k = 200;
for k = 1:64
    subplot(8, 8, k);
%     trisurf(faces, verts(:, 1), verts(:, 2), verts(:, 3), W(:, k), 'EdgeColor', 'none'); view(2); axis image off; shading interp;
    trisurf(faces, verts(:, 1), verts(:, 2), verts(:, 3), V(:, k), 'EdgeColor', 'none'); view(2); axis image off; shading interp; colormap viridis;
end



end