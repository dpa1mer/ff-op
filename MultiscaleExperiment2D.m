function [meanEdgeLengths, evs] = MultiscaleExperiment2D(meshfile, levels)

[verts{1}, faces{1}] = load_mesh(meshfile);
meshData{1} = ProcessMesh2D(verts{1}, faces{1});
verts{1} = meshData{1}.verts;
faces{1} = meshData{1}.faces;

%% Build mesh at each level
for level = 2:levels
    [verts{level}, faces{level}, P{level-1}, I{level-1}] = loop(verts{level-1}, faces{level-1});
    meshData{level} = ProcessMesh2D(verts{level}, faces{level});
end

%% Compute frame field at finest level and average up
z4{levels} = MBO2D(meshData{levels}, true);
for level = (levels-1):-1:1
    z4{level} = accumarray(I{level}, z4{level + 1}, [meshData{level}.nf, 1]);
    z4{level} = z4{level} ./ abs(z4{level});
end

for level = 1:levels
    meanEdgeLengths(level) = mean(meshData{level}.edgeLengths);
    fprintf('Level %d: mean edge length %d\n', level, meanEdgeLengths(level));
    [~, Tij] = Frame2Tensor2D(meshData{level}, z4{level});
    [Op, M] = PhaseField2D(meshData{level}, Tij, true);
    [V, D] = eigs(Op + 1e-6 * M, M, 64, 'smallestabs');
    figure('Name', sprintf('Level %d', level));
    trisurf(meshData{level}.faces, meshData{level}.verts(:, 1), meshData{level}.verts(:, 2), meshData{level}.verts(:, 3), V(:, 64) * sign(V(100, 64)), 'EdgeColor', 'none', 'FaceColor', 'flat'); view(2); axis image off; shading interp; colormap viridis;
    evs(:, level) = diag(D);
end

figure('Name', 'Eigenvalue Convergence');
bar(evs(:, 1:levels-1) - evs(:, levels));
legend('Location', 'northwest');

end
