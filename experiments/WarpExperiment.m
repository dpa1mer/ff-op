function evs = WarpExperiment(ellipticity)

if nargin < 1
    ellipticity = 1e-2;
end

% fd=@(p) dunion(drectangle(p,-1,1,-1,0),drectangle(p,-0.25,1,-0.5,1));
% figure; [p,faces]=distmesh2d(fd,@huniform,0.02,[-1,-1;1,1],[-1 -1; -1 0; -0.25 0; -0.25 1; 1 1; 1 -1]);%; 0.5 -0.5; 0.5 -1]);
% verts = [p, zeros(size(p, 1), 1)];
[verts, faces] = load_mesh('../../models/polyrect2.obj');
figure; trisurf(faces, verts(:, 1), verts(:, 2), verts(:, 3), sin(8*pi*verts(:, 1)) .* sin(8*pi*verts(:, 2)), 'EdgeColor', 'none'); view(2); axis image off; shading interp; colormap([0 0 0; 1 0 1]);
nv = size(verts, 1);
domain{1} = ProcessMesh2D(verts, faces);
Tij{1} = repmat(diag([ellipticity, 1, 1, ellipticity]), 1, 1, nv);

%% Warp the vertex positions
z = verts(:, 1) + 1i*verts(:, 2);
f = @(z) (z + 2 - 2i).^2;
df = @(z) 2 .* (z + 2 - 2i);
fz = f(z);
verts2 = [real(fz), imag(fz), verts(:, 3)];
figure; trisurf(faces, verts2(:, 1), verts2(:, 2), verts2(:, 3), sin(8*pi*verts(:, 1)) .* sin(8*pi*verts(:, 2)), 'EdgeColor', 'none'); view(2); axis image off; shading interp; colormap([0 0 0; 1 0 1]);
domain{2} = ProcessMesh2D(verts2, faces);

%% Compute the map differential frame
u = df(z);
v = 1i * u;
uv = [reshape([real(u), imag(u)].', 2, 1, []), reshape([real(v), imag(v)].', 2, 1, [])];

xy = uv;%multitransp(batchop('pinv', uv, 2));
x2 = reshape(reshape(xy(:,1,:), 2, 1, nv) .* reshape(xy(:,1,:), 1, 2, nv), 4, 1, nv);
x4 = x2 .* reshape(x2, 1, 4, nv);
y2 = reshape(reshape(xy(:,2,:), 2, 1, nv) .* reshape(xy(:,2,:), 1, 2, nv), 4, 1, nv);
y4 = y2 .* reshape(y2, 1, 4, nv);
Tij{2} = reshape(max(batchop('eig', x4 + y4), [], 1), 1, 1, []) .* eye(4) - (1 - ellipticity) * (x4 + y4);

% figure; FancyQuiver([verts2;verts2], [[squeeze(xy(:, 1, :)).'; squeeze(xy(:, 2, :)).'] zeros(2*nv, 1)], viridis, 0); view(2); axis image off;

%% Compute and plot eigenfunctions on each domain
for j = 1:2
    [Op, M] = PhaseField2D(domain{j}, Tij{j}, true);
    [V, D] = eigs(Op + 1e-6 * M, M, 65, 'smallestabs');
    evs(:, j) = diag(D);
    figure; tlt = tiledlayout(8, 8);
    tlt.TileSpacing = 'none';
    [cmin, cmax] = bounds(V, 'all');
    for k = 1:64
        nexttile(k);
        trisurf(faces, verts2(:, 1), verts2(:, 2), verts2(:, 3), V(:, k + 1) .* sign(V(100, k + 1)), 'EdgeColor', 'none');
        view(2); axis image off; shading interp; colormap viridis; %caxis([cmin, cmax]);
    end
end

end