addpath diffusioncurves;

[V,F,b,bc] = flower_closed();
us = DiffusionCurvesClosed(V, F, b, bc, 'biharmonic');

hold on;
for i = 1:numel(V)
    tsurf(F{i}, V{i}, 'FaceVertexCData',us{i});
end
hold off;
axis equal;
shading interp;
