function [crossFaceTensor, crossVtxTensor] = Frame2Tensor2D(meshData, z4, ellipticity)

if nargin < 3
    ellipticity = 0.01;
end

nf = meshData.nf;
nv = meshData.nv;

z = (z4.^(1/4)) .* [1 1i -1 -1i];

% Form cross field tensor
u = [real(z(:, 1)), imag(z(:, 1))].';
v = [real(z(:, 2)), imag(z(:, 2))].';
u2 = reshape(reshape(u, 2, 1, nf) .* reshape(u, 1, 2, nf), 4, 1, nf);
u4 = u2 .* reshape(u2, 1, 4, nf);
v2 = reshape(reshape(v, 2, 1, nf) .* reshape(v, 1, 2, nf), 4, 1, nf);
v4 = v2 .* reshape(v2, 1, 4, nf);
crossFaceTensor = eye(4) - (1 - ellipticity) * (u4 + v4);

% Average T to the vertices
areaWeightedT = repmat(reshape(meshData.areas, 1, 1, nf) .* crossFaceTensor, 1, 1, 1, 3) / 3;
MT_I = repmat([1;2;3;4], 1, 4, nf * 3);
MT_J = repmat([1 2 3 4], 4, 1, nf * 3);
MT_K = repmat(reshape(meshData.faces(:), 1, 1, nf * 3), 4, 4, 1);
MTij = accumarray([MT_I(:) MT_J(:) MT_K(:)], areaWeightedT(:), [4 4 nv]);
crossVtxTensor = MTij ./ reshape(full(diag(meshData.star0lump)), 1, 1, nv);

end