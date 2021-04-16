function us = DiffusionCurvesClosed(V, F, b, bc, operator)
%DIFFUSIONCURVESCLOSED On a mesh V,F solve divvusion curve problem with boundary
%vertices b, boundary colors bc, and operator operator.
%Each of the arguments is a cell of things that should be composed

if nargin<5
    operator = 'biharmonic';
end

us = {};

%Find domain for each th and mesh and solve
for i = 1:numel(V)
    cV = V{i};
    cF = F{i};
    cb = b{i};
    cbc = bc{i};
    
    %Solve problem
    switch operator
        case 'biharmonic'
            L = -cotmatrix(cV,cF);
            M = massmatrix(cV,cF,'voronoi');
            Q = L' * (M\L);
            u = min_quad_with_fixed(0.5*Q, zeros(size(cV,1),1), cb, cbc);
    end
    us{end+1} = u;
end

end

