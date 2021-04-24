function [Op, star0lump] = FFOp3D(meshData, q, ellipticity, neumann)

if nargin < 4
    neumann = true;
end

if nargin < 3
    ellipticity = 0.01;
end

if nargin < 2
    %% Compute frame field by MBO

    q = MBO(meshData, RayMBO);
    octa = OctaMBO;
    q = octa.proj(q);
end

%% Load mesh.
assert(~meshData.isGrid);

verts = meshData.verts;
tets = meshData.tets;
nt = meshData.nt;
nv = meshData.nv;
bdryIdx = meshData.bdryIdx;
nb = length(bdryIdx);

%% Convert frames to tensors
frameTensors = reshape(Monomial2Tensor(Sph024ToMonomial(Octa2Odeco(q))), 9, 9, nv);
Tij = eye(9) - (1 - ellipticity) * (frameTensors ./ 24);

%% Construct phase field operators
volumes = TetVolumes(verts, tets);

tetOppFaces = reshape(tets(:, flipud(nchoosek(1:4, 3))), nt, 4, 3);
oppFaceNormals = cross(verts(tetOppFaces(:, :, 2), :) - verts(tetOppFaces(:, :, 1), :), ...
                       verts(tetOppFaces(:, :, 3), :) - verts(tetOppFaces(:, :, 1), :), 2);
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

if neumann
    %% Impose 0-Neumann Boundary Conditions (weakly)
    bdryNormals = reshape(meshData.bdryNormals', 3, 1, nb);
    nn = batchop('mult', bdryNormals, bdryNormals, 'N', 'T');
    nnn = batchop('mult', bdryNormals, reshape(nn, 1, 9, nb));
    multN = [multitransp(bdryNormals) zeros(1, 6, nb);
             zeros(1, 3, nb) multitransp(bdryNormals) zeros(1, 3, nb);
             zeros(1, 6, nb) multitransp(bdryNormals)];
    multNt = multitransp(reshape(multN, 9, 3, nb));
    Bij = [multN - nnn; multNt - nnn];
    BT = batchop('mult', Bij, Tij(:,:,bdryIdx));
    BTB = batchop('mult', BT, Bij, 'N', 'T');
    BTBi = batchop('pinv', BTB, 4);
    Tij(:,:,bdryIdx) = Tij(:,:,bdryIdx) - batchop('mult', batchop('mult', BT, BTBi, 'T', 'N'), BT, 'N', 'N');
end

%% Form frame field principal symbol matrix
vertBaseIdx = reshape(9 * (0:nv - 1), 1, 1, nv);
T_I = repmat(vertBaseIdx + (1:9).', 1, 9, 1, 1);
T_J = repmat(vertBaseIdx + (1:9), 9, 1, 1, 1);
T = sparse(T_I(:), T_J(:), Tij(:), 9 * nv, 9 * nv);
TM = T / M9;

if ~neumann
    %% Only take interior nodes
    intBlkIdx = (9 * meshData.intIdx + (-8:0)).';
    D = D(:, intBlkIdx);
    TM = TM(intBlkIdx, intBlkIdx);
end

%% Build final operator
DVG = D' * V * G;
Op = DVG' * TM * DVG;

end