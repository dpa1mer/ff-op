function [us,zs] = DiffusionCurvesClosed(V, F, b, bc, operator, rotfield)
%DIFFUSIONCURVESCLOSED On a mesh V,F solve diffusion curve problem with boundary
%vertices b, boundary colors bc, and operator operator.
%Each of the arguments is a cell of things that should be composed

if nargin<5
    operator = 'biharmonic';
end
if nargin<6
    rotfield = false;
end

us = {};
zs = {};

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
            
            u = nan(size(cV,1),3);
            for j=1:3
                u(:,j) = quadprog(Q, zeros(size(cV,1),1), [], [], ...
                    sparse(1:numel(cb),cb,ones(numel(cb),1),numel(cb),size(cV,1)), cbc(:,j), ...
                    min(cbc(:,j))*ones(size(cV,1),1), ...
                    max(cbc(:,j))*ones(size(cV,1),1));
            end
            u = min_quad_with_fixed(0.5*Q, zeros(size(cV,1),1), cb, cbc);
        case 'odeco'
            %Create the Palmer data structure
            mesh = ProcessMesh2D([cV zeros(size(cV,1),1)], cF);
            
            %Create operator
            z = MBO2D(mesh, false);
            if rotfield
                z = z * exp(1i*pi);
            end
            [~, Tij] = Frame2Tensor2D(mesh, z, 0.01);
            [Op, M] = PhaseField2D(mesh, Tij, false);
            
            u = nan(size(Op,1),3);
            for j=1:3
                u(:,j) = quadprog(Op, zeros(size(cV,1),1), [], [], ...
                    sparse(1:numel(cb),cb,ones(numel(cb),1),numel(cb),size(cV,1)), cbc(:,j), ...
                    min(cbc(:,j))*ones(size(cV,1),1), ...
                    max(cbc(:,j))*ones(size(cV,1),1));
            end
            
            zs{end+1} = z;
    end
    u(cb,:) = cbc;
    us{end+1} = u;
end

end

