function weights = BOWExperiment(meshData, nh)

verts = meshData.verts;
faces = meshData.faces;
nv = meshData.nv;

% Use neumann BCs
[~, Tij] = Frame2Tensor2D(meshData, MBO2D(meshData, true), 1e-2);
[Op, ~] = PhaseField2D(meshData, Tij, true);

% Select random handles
handleIdx = randi(nv, [nh 1]);

cvx_begin
    cvx_precision high

    variable weights(nv, nh)
    
    qf = quad_form(weights(:, 1), Op);
    for i = 2:nh
        qf = qf + quad_form(weights(:, i), Op);
        i
    end
    
    minimize qf
    subject to
        weights(handleIdx, :) == eye(nh);
        sum(weights, 2) == 1
        0 <= weights <= 1
cvx_end

figure;
for k = 1:nh
    subplot(8, 8, k);
    trisurf(faces, verts(:, 1), verts(:, 2), verts(:, 3), weights(:, k), 'EdgeColor', 'none'); view(2); axis image off; shading interp; colormap viridis;
    hold on; scatter3(verts(handleIdx, 1), verts(handleIdx, 2), verts(handleIdx, 3), 200, 'r.');
end

end