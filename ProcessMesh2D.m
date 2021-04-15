function meshData = ProcessMesh2D(verts, faces)

faces = flip_ears(verts, faces);

meshData.verts = verts; meshData.faces = faces;
meshData.faceOrientedEdges = reshape(faces(:, [1 2; 2 3; 3 1]), [], 2);
meshData.faceOrientation = reshape(sign(meshData.faceOrientedEdges(:, 2) - meshData.faceOrientedEdges(:, 1)), [], 3);
faceEdges = sort(meshData.faceOrientedEdges, 2);
[meshData.edges, ~, meshData.face2edge] = unique(faceEdges, 'rows');
meshData.face2edge = reshape(meshData.face2edge, [], 3);
ne = size(meshData.edges, 1); meshData.ne = ne;
nv = size(verts, 1); meshData.nv = nv;
nf = size(faces, 1); meshData.nf = nf;

meshData.edgeLengths = vecnorm(meshData.verts(meshData.edges(:, 2), :) - meshData.verts(meshData.edges(:, 1), :), 2, 2);
meshData.faceCenters = (1/3) * squeeze(sum(reshape(verts(faces, :), nf, 3, 3), 2));

%% Construct Hodge Stars
vertOppEdges = meshData.face2edge(:, [2 3 1]);
eij = verts(faces(:, [3 1 2]), :) - verts(faces, :);
eik = verts(faces(:, [2 3 1]), :) - verts(faces, :);
meshData.faceVertAngles = reshape(acos((dot(eij, eik, 2) ./ vecnorm(eij, 2, 2)) ./ vecnorm(eik, 2, 2)), nf, 3);
oppositeCotans = cot(meshData.faceVertAngles);

Lij = oppositeCotans / 2;
L = sparse(meshData.edges(vertOppEdges, 1), meshData.edges(vertOppEdges, 2), Lij, nv, nv);
meshData.star1 = spdiags(nonzeros(L'), 0, ne, ne);
meshData.star1dual = spdiags(diag(meshData.star1).^(-1), 0, ne, ne);

meshData.areas = 0.5 * vecnorm(cross(verts(faces(:, 1), :) - verts(faces(:, 2), :), verts(faces(:, 3), :) - verts(faces(:, 2), :)), 2, 2);
meshData.star2 = spdiags(meshData.areas, 0, nf, nf);

meshData.star0 = massmatrix(verts, faces, 'full');%sparse(faces, faces, repmat(areas / 3, 1, 3), nv, nv);
meshData.star0lump = massmatrix(verts, faces, 'barycentric');

%% Construct Laplacians
meshData.d0 = sparse(repmat((1:ne).', 1, 2), meshData.edges, repmat([-1 1], ne, 1), ne, nv);
meshData.d1 = sparse(repmat((1:nf).', 1, 3), meshData.face2edge, meshData.faceOrientation, nf, ne);
meshData.L = meshData.d0.' * meshData.star1 * meshData.d0;
meshData.Ldual = meshData.d1 * meshData.star1dual * meshData.d1.';

%% Construct boundary matrices
tri = triangulation(faces, verts);
meshData.bdryEdges = freeBoundary(tri);
meshData.bdryIdx = unique(meshData.bdryEdges);
meshData.intIdx = setdiff(1:meshData.nv, meshData.bdryIdx);
nb = size(meshData.bdryEdges, 1); meshData.nb = nb;

[bdryForward, bdryEdgeIdxF] = ismember(meshData.bdryEdges, meshData.edges, 'rows');
[bdryBackward, bdryEdgeIdxB] = ismember(fliplr(meshData.bdryEdges), meshData.edges, 'rows');
meshData.bdryEdgeIdx = bdryEdgeIdxF + bdryEdgeIdxB;

bdryEdgeVecs = verts(meshData.bdryEdges(:, 2), :) - verts(meshData.bdryEdges(:, 1), :);
meshData.edgeNormals = [-bdryEdgeVecs(:, 2), bdryEdgeVecs(:, 1), bdryEdgeVecs(:, 3)];
meshData.edgeNormals = meshData.edgeNormals ./ vecnorm(meshData.edgeNormals, 2, 2);
meshData.bdryEdgeCenters = 0.5 * (verts(meshData.bdryEdges(:, 2), :) + verts(meshData.bdryEdges(:, 1), :));

vtxNormals = accumarray([repmat(meshData.bdryEdges(:, 1), 3, 1) repelem((1:3).', nb)], meshData.edgeNormals(:), [nv 3]) ...
           + accumarray([repmat(meshData.bdryEdges(:, 2), 3, 1) repelem((1:3).', nb)], meshData.edgeNormals(:), [nv 3]);
vtxNormals = vtxNormals(meshData.bdryIdx, :);
meshData.vtxNormals = vtxNormals ./ vecnorm(vtxNormals, 2, 2);

end