function w = QuasiEig(meshData, lambda)

verts = meshData.verts;
faces = meshData.faces;
nv = meshData.nv;

[~, Tij] = Frame2Tensor2D(meshData, MBO2D(meshData, true), 1e-2);
[Op, M, Bnn] = FFOp2D(meshData, Tij, true);
BigOp = [(Op - lambda * M)' * (Op - lambda * M), Bnn'; Bnn zeros(size(Bnn, 1))];
w = BigOp \ [zeros(nv, 1); -ones(size(Bnn, 1), 1)];
w = w(1:nv);
figure; trisurf(faces, verts(:, 1), verts(:, 2), verts(:, 3), w, 'EdgeColor', 'none');
view(2); axis image off; shading interp; colormap viridis;

end