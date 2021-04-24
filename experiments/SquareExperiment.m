function [evsAnalytic, evsEmpirical, meanEdgeLengths] = SquareExperiment

%% Compute spectrum analytically
[Kx, Ky] = meshgrid(0:200, 0:200);
Kx = Kx .* (pi/2); Ky = Ky .* (pi/2);
evsAnalytic = sort(1e-2 .* (Kx(:).^4 + Ky(:).^4) + 2 .* Kx(:).^2 .* Ky(:).^2);
figure; plot(evsAnalytic(1:100)); hold on;

%% Compute spectrum empirically
meshFiles = ["meshes/square0.03.obj", "meshes/square0.025.obj", ...
             "meshes/square0.02.obj", "meshes/square0.015.obj"];

i = 1;
for meshFile = meshFiles
    [verts, faces] = load_mesh(meshFile);
    square = ProcessMesh2D(verts, faces);
    meanEdgeLengths(i) = mean(square.edgeLengths);
    [~, Tij] = Frame2Tensor2D(square, MBO2D(square, true), 1e-2);
    [Op, M] = FFOp2D(square, Tij, true);
    evsEmpirical(:, i) = eigs(Op + 1e-6 * M, M, 100, 'smallestabs');
    i = i + 1;
end
plot(evsEmpirical);


end