function EllipticityExperiment2D(ellipticity, diffusionTime)

[verts, faces] = load_mesh('../../models/disk2.obj');
disk = ProcessMesh2D(verts, faces);
nv = disk.nv;
[~, Tij] = Frame2Tensor2D(disk, MBO2D(disk, true), ellipticity);
[Op, M] = PhaseField2D(disk, Tij, true);

%% Green's Functions
vIdx = [3000, 4000, 5000, 6000, 9000, 11000, 13000, 15000];
delta = zeros(nv, 1);
delta(vIdx) = 1;

green = (M + diffusionTime * Op) \ (M * delta);
figure; trisurf(faces, verts(:, 1), verts(:, 2), verts(:, 3), green, 'EdgeColor', 'none');
view(2); axis image off; shading interp; colormap viridis;
hold on; scatter3(verts(vIdx, 1), verts(vIdx, 2), verts(vIdx, 3), 200, 'r.');

end