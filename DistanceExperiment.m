function DistanceExperiment(meshData, ellipticity, bcNeumann, diffusionTime)

if nargin < 3
    bcNeumann = false;
end

if nargin < 2
    ellipticity = 1e-2;
end

nv = meshData.nv;
verts = meshData.verts;
faces = meshData.faces;

[~, Tij] = Frame2Tensor2D(meshData, MBO2D(meshData, true), ellipticity);

% Use 0-Neumann BCs
[Op, M] = PhaseField2D(meshData, Tij, bcNeumann);

%% Approximate distance using eigenbasis
[W, D] = eigs(Op + 1e-6 * M, M, 200, 'smallestabs', 'IsSymmetricDefinite', true);
% Get rid of kernel
if bcNeumann
    W = W(:, 2:end); D = D(2:end, 2:end);
else
    W = W(:, 4:end); D = D(4:end, 4:end);
end

vIdx = 2567;%4550;%21084;%randi(length(verts), 1);
dist = sqrt(sum((W - W(vIdx, :)).^2 ./ (diag(D).').^2, 2));
figure; trisurf(faces, verts(:, 1), verts(:, 2), verts(:, 3), dist, 'EdgeColor', 'none');
view(2); axis image off; shading interp; colormap viridis(20);
hold on; scatter3(verts(vIdx, 1), verts(vIdx, 2), verts(vIdx, 3), 200, 'r.');

% %% Green's Functions
% vIdx = randi(length(verts), [10, 1]);
% delta = zeros(nv, 1);
% delta(vIdx) = 1;
% 
% green = (M + diffusionTime * Op) \ delta;
% figure; trisurf(faces, verts(:, 1), verts(:, 2), verts(:, 3), green, 'EdgeColor', 'none');
% view(2); axis image off; shading interp; colormap viridis;
% hold on; scatter3(verts(vIdx, 1), verts(vIdx, 2), verts(vIdx, 3), 200, 'r.');
% 
% greenProj = W * (W' * M * green);
% figure; trisurf(faces, verts(:, 1), verts(:, 2), verts(:, 3), greenProj, 'EdgeColor', 'none');
% view(2); axis image off; shading interp; colormap viridis;
% hold on; scatter3(verts(vIdx, 1), verts(vIdx, 2), verts(vIdx, 3), 200, 'r.');

end