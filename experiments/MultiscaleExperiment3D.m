function [meanEdgeLengths, evs] = MultiscaleExperiment3D(colormap)

% meshfiles = ["meshes/spheres/sphere_r0.15.mesh";
%              "meshes/spheres/sphere_r0.14.mesh";
%              "meshes/spheres/sphere_r0.13.mesh";
%              "meshes/spheres/sphere_r0.12.mesh";
%              "meshes/spheres/sphere_r0.11.mesh";
%              "meshes/spheres/sphere_r0.10.mesh";
%              "meshes/spheres/sphere_r0.09.mesh";
%              "meshes/spheres/sphere_r0.08.mesh";
%              "meshes/spheres/sphere_r0.07.mesh";
%              "meshes/spheres/sphere_r0.06.mesh";
%              "meshes/spheres/sphere_r0.05.mesh";
%              "meshes/spheres/sphere_r0.04.mesh"];
% %              "meshes/spheres/sphere_r0.03.mesh";
% %              "meshes/spheres/sphere_r0.02.mesh"];

meshfiles = ["../meshes/dbommes-meshes/teddy_80k.mesh";
             "../meshes/dbommes-meshes/teddy_200k.mesh";
             "../meshes/dbommes-meshes/teddy_370k.mesh"];

levels = length(meshfiles);

%% Build mesh at each level
for level = 1:levels
    meshData{level} = ImportMesh(meshfiles{level});
end

%% Compute frame field at finest level and average up
q{levels} = MBO(meshData{levels}, RayMBO, [], 50, 3);
octa = OctaMBO;
q{levels} = octa.proj(q{levels});
qT = q{levels}';

for level = (levels-1):-1:1
    [tetIx, bary] = pointLocation(meshData{levels}.tetra, meshData{level}.verts);
    insideIx = isfinite(tetIx);
    q{level} = zeros(9, meshData{level}.nv);
    q{level}(:, insideIx) = squeeze(sum(bary(insideIx, :) .* reshape(qT(meshData{levels}.tets(tetIx(insideIx), :), :), [], 4, 9), 2))';
    q{level}(:, insideIx) = octa.proj(q{level}(:, insideIx));
    q{level} = MBO(meshData{level}, RayMBO, q{level}, 1, 3);
    q{level} = octa.proj(q{level});
end

for level = 1:levels
        meanEdgeLengths(level) = mean(vecnorm(meshData{level}.verts(meshData{level}.edges(:, 2), :) - ...
                                              meshData{level}.verts(meshData{level}.edges(:, 1), :), 2, 2));
        fprintf('Level %d: mean edge length %d\n', level, meanEdgeLengths(level));
        [Op{level}, M{level}] = FFOp3D(meshData{level}, q{level});
end

%% Solve eigenvalue problem
for level = 1:levels
    [V, D] = eigs(Op{level} + 1e-6 * M{level}, M{level}, 64, 'smallestabs');
    [~, maxXIdx] = max(meshData{level}.verts(:, 1)); % Consistent sign: make EFs positive here
    fig = figure('Name', sprintf('Level %d : mean edge length %d', level, meanEdgeLengths(level)), 'NumberTitle', 'off');
    FancyScatter(meshData{level}.verts, V(:, 32) * sign(V(maxXIdx, 32)), colormap, 0, 2); view(30, 70); axis image vis3d off;
    fig.WindowState = 'maximized';
    evs(:, level) = diag(D);
end

VisualizeFrameField(meshData{levels}, q{levels});

figure('Name', 'Eigenvalue Convergence');
plot(evs);
legend('Location', 'northwest');

end
