function V = PhaseField3DGrid(grid, q)

s = sqrt(2);

tens = Monomial2Mandel(Sph024ToMonomial(Octa2Odeco(q)));
tens = eye(6) - tens;
gridHess = [grid.Dxx; grid.Dyy; grid.Dzz; s * grid.Dyz; s * grid.Dxz; s * grid.Dxy];
ii = repmat(grid.nv * (0:5).', 1, 6, grid.nv) + reshape(1:grid.nv, 1, 1, grid.nv);
jj = repmat(grid.nv * (0:5), 6, 1, grid.nv) + reshape(1:grid.nv, 1, 1, grid.nv);
Orthotropic = gridHess' * sparse(ii(:), jj(:), tens(:), 6 * grid.nv, 6 * grid.nv) * gridHess;
[V, D] = eigs(Orthotropic + 1e-6 * speye(grid.nv), 200, 'smallestabs');

fig = figure; colormap inferno;
tiledlayout(5, 5);
for j = 1:25
    ix = V(:, j).^2 > mean(V(:, j).^2);
    ax(j) = nexttile(j); scatter3(grid.verts(ix, 1), grid.verts(ix, 2), grid.verts(ix, 3), 1, V(ix, j), '.');
    view(3);
    axis image vis3d off;
end
% linkaxes(ax, 'xyz');
Link = linkprop(ax, {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
setappdata(fig, 'TheLink', Link);

end