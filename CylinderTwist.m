function cylinder = CylinderTwist(cylinder)

bdryIdx = cylinder.bdryIdx;
bdryNormals = cylinder.bdryNormals;
verts = cylinder.verts;
nb = length(bdryIdx);

% Extract cap and shell
capIdx = find(vecnorm(cross(bdryNormals, repmat([0 0 1], [nb, 1]), 2), 2, 2) < 0.9);
shellIdx = setdiff((1:nb)', capIdx);

% Make cap unconstrained
bdryIdx = bdryIdx(shellIdx);
bdryNormals = bdryNormals(shellIdx, :);

% Add constraints in a shell
constrainedIdx = find((vecnorm(verts(:, 1:2), 2, 2) > 50) & (vecnorm(verts(:, 1:2), 2, 2) < 60));
constrainedIdx = setdiff(constrainedIdx, bdryIdx);
constrainedVerts = verts(constrainedIdx, :);
bdryIdx = [bdryIdx; constrainedIdx];

% Compute constraint direction
constraintNormals = [(1 ./ vecnorm(constrainedVerts(:, 1:2), 2, 2)) .* [-constrainedVerts(:, 2), constrainedVerts(:, 1)], 2*ones(length(constrainedIdx), 1)];
constraintNormals = constraintNormals ./ vecnorm(constraintNormals, 2, 2);
bdryNormals = [bdryNormals; constraintNormals];

cylinder.bdryIdx = bdryIdx;
cylinder.bdryNormals = bdryNormals;
cylinder.intIdx = setdiff((1:cylinder.nv)', cylinder.bdryIdx);

end