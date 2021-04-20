function [meanEdgeLengths, evs] = MultiscaleExperiment3D()

meshfiles = ["meshes/spheres/sphere_r0.15.mesh";
             "meshes/spheres/sphere_r0.14.mesh";
             "meshes/spheres/sphere_r0.13.mesh";
             "meshes/spheres/sphere_r0.12.mesh";
             "meshes/spheres/sphere_r0.11.mesh";
             "meshes/spheres/sphere_r0.10.mesh";
             "meshes/spheres/sphere_r0.09.mesh";
             "meshes/spheres/sphere_r0.08.mesh";
             "meshes/spheres/sphere_r0.07.mesh";
             "meshes/spheres/sphere_r0.06.mesh";
             "meshes/spheres/sphere_r0.05.mesh"];

levels = length(meshfiles);

%% Build mesh at each level
for level = 1:levels
    meshData{level} = ImportMesh(meshfiles{level});
end

%% Compute frame field at finest level and average up
q{levels} = MBO(meshData{levels}, RayMBO);
octa = OctaMBO;
q{levels} = octa.proj(q{levels});
qT = q{levels}';

for level = (levels-1):-1:1
    [tetIx, bary] = pointLocation(meshData{levels}.tetra, meshData{level}.verts);
    insideIx = isfinite(tetIx);
    q{level} = zeros(9, meshData{level}.nv);
    q{level}(:, insideIx) = squeeze(sum(bary(insideIx, :) .* reshape(qT(meshData{levels}.tets(tetIx(insideIx), :), :), [], 4, 9), 2))';
    q{level}(:, insideIx) = octa.proj(q{level}(:, insideIx));
    q{level} = MBO(meshData{level}, RayMBO, q{level});
    q{level} = octa.proj(q{level});
end

for level = 1:levels
        meanEdgeLengths(level) = mean(vecnorm(meshData{level}.verts(meshData{level}.edges(:, 2), :) - ...
                                              meshData{level}.verts(meshData{level}.edges(:, 1), :), 2, 2));
        fprintf('Level %d: mean edge length %d\n', level, meanEdgeLengths(level));
        [Op{level}, M{level}] = PhaseField3D(meshData{level}, q{level});
end

%% Solve eigenvalue problem
for level = 1:levels
    [V, D] = eigs(Op{level} + 1e-6 * M{level}, M{level}, 64, 'smallestabs');
    figure('Name', sprintf('Level %d : mean edge length %d', level, meanEdgeLengths(level)), 'NumberTitle', 'off');
    plot_3df(meshData{level}.verts, meshData{level}.tets, V(:, 64), 'PlotType', 'CleanSlice');
    evs(:, level) = diag(D);
end

figure('Name', 'Eigenvalue Convergence');
bar(evs(:, 1:levels-1) - evs(:, levels));
legend('Location', 'northwest');

end
