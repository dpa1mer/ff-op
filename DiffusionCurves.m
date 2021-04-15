function us = DiffusionCurves(V, F, b, bc)
%DIFFUSION_CURVES On a mesh V,F solve divvusion curve problem with boundary
%vertices b, boundary colors bc.
%Each of the arguments is a cell of things that should be composed

us = {};

%Find domain for each th and mesh and solve
hold on;
for i = 1:numel(V)
    cV = V{i};
    cF = F{i};
    cb = b{i};
    cbc = bc{i};
    
    %Solve problem
    %Solve biharmonic interpolation problem
    L = -cotmatrix(cV,cF);
    M = massmatrix(cV,cF,'voronoi');
    Q = L' * (M\L);
    keyboard
    u = min_quad_with_fixed(0.5*Q, zeros(size(cV,1),1), cb, cbc);
    us{end+1} = u;
    
    tsurf(cF, cV, 'FaceVertexCData',u);
end
hold off;
axis equal;
shading interp;

end

