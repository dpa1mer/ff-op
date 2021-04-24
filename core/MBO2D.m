function z4 = MBO2D(meshData, dual, tauMult)

if nargin < 2
    dual = false;
end

if nargin < 3
    tauMult = 0.5;
end

nv = meshData.nv;
nf = meshData.nf;

%% Set up optimization and boundary conditions
if dual
    n = nf;
    zb = (meshData.edgeNormals(:, 1) + 1i.*meshData.edgeNormals(:, 2)).^4;
    B = abs(meshData.d1.');
    B = B(meshData.bdryEdgeIdx, :);
    L = meshData.Ldual;
    M = meshData.star2;
else
    n = nv;
    zb = (meshData.vtxNormals(:, 1) + 1i.*meshData.vtxNormals(:, 2)).^4;
    B = sparse(1:meshData.nb, meshData.bdryIdx, 1, meshData.nb, meshData.nv);
    L = meshData.L;
    M = meshData.star0;
end

%% Compute cross field by MBO

lambda = eigs(L, M, 2, 'smallestabs');
tau = tauMult/lambda(2);

A = [M + tau * L, B.'; B zeros(size(B, 1))];
cholA = decomposition(A, 'ldl');
z4 = zeros(n, 1);
for j = 1:10
    z4 = cholA \ [M * z4; zb];
    z4 = z4(1:n);
    crossNorm = abs(z4);
    z4 = z4 ./ crossNorm;
end

end