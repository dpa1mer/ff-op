function EllipticityExperiment2D(meshData, nDeltas, diffusionTime)

if nargin < 3
    diffusionTime = 1e-5;
end

if nargin < 2
    nDeltas = 16;
end

verts = meshData.verts;
faces = meshData.faces;
nv = meshData.nv;
z4 = MBO2D(meshData);
% z4 = (meshData.verts(:, 1) + 1i * meshData.verts(:, 2));
% z4 = (z4 ./ abs(z4)).^4 .* exp(2*pi*1i*(1 - abs(z4)));
figure('Name', 'Frame Field', 'WindowState', 'maximized'); VisualizeFrameField2D(meshData, z4);

%% Generate good samples
meanEdgeLength = mean(meshData.edgeLengths);
minDist = 0;
while minDist < 15 * meanEdgeLength
    vtxWeights = diag(meshData.star0lump);
    vIdx = meshData.intIdx(randsample(length(meshData.intIdx), nDeltas, true, vtxWeights(meshData.intIdx)));
    gram = meshData.verts(vIdx, :) * meshData.verts(vIdx, :)';
    dist = sqrt(diag(gram) + diag(gram)' - 2 * gram);
    dist = dist + diag(inf(nDeltas, 1)); % Ignore diagonal
    minDist = min(dist, [], 'all');
end

delta = zeros(nv, 1);
delta(vIdx) = 1;

%% Compute Impulse Responses
for i = 0:-1:-3
    ellipticity = 5^i;
    [~, Tij] = Frame2Tensor2D(meshData, z4, ellipticity);
    [Op, M] = PhaseField2D(meshData, Tij, false); % Natural boundary conditions
    green = (M + diffusionTime * Op) \ (M * delta);
    figure('Name', sprintf('Ellipticity 5^%d', i), 'WindowState', 'maximized');
    trisurf(faces, verts(:, 1), verts(:, 2), verts(:, 3), green, 'EdgeColor', 'none');
    view(2); axis image off; shading interp; colormap(cbrewer('YlGnBu', 500));
    hold on; scatter3(verts(vIdx, 1), verts(vIdx, 2), verts(vIdx, 3), 200, 'c.');
end

end