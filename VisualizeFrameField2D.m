function VisualizeFrameField2D(meshData, z4)

verts = meshData.verts;
faces = meshData.faces;
tr = triangulation(faces, verts(:, 1:2));

z = (z4.^(1/4)) .* [1 1i -1 -1i];

np = 2000;
points = verts(meshData.intIdx(randi(length(meshData.intIdx), np, 1)), 1:2); % start one path at each interior vertex
np = size(points,1); % number of points moving around
vel = randn(np, 2); % initial direction

% stepsize for gradient descent
stepsize = 0.5 * (mean(meshData.edgeLengths) - std(meshData.edgeLengths));

[lb, ub] = bounds(meshData.verts, 1);
curveLength = 0.25 * max(ub - lb);

% number of steps to flow
nsteps = round(curveLength / stepsize);

paths = nan(np, 2, nsteps);
paths(:, :, 1) = points;

% Trace integral curves of frame field
for i=1:nsteps
    % Interpolate frame field
    vel = interpVelocity(points, vel);
    
    % Normalize velocity
    vel = vel ./ vecnorm(vel, 2, 2);
    
    % Update points
    points = points + stepsize * vel; 
    
    % Add points to paths
    paths(:, :, i) = points;
end

figure; plot(squeeze(paths(:, 1, :))', squeeze(paths(:, 2, :))', 'LineWidth', 2);
axis image off;
colororder(cbrewer('Paired', 12));

% Returns interpolation of the frame field component best matching old velocity
function newVel = interpVelocity(pts, oldVel)
    % Use complex numbers
    oldVel = oldVel(:, 1) + 1i * oldVel(:, 2);
    
    % Get barycentric coords
    badpts = any(isnan(pts), 2);
    fIdx = ones(size(pts, 1), 1);
    bary = zeros(size(pts, 1), 3);
    [fIdx(~badpts), bary(~badpts, :)] = pointLocation(tr, pts(~badpts, :));
    badpts = badpts | isnan(fIdx);
    fIdx(badpts) = 1; % placeholder index
    
    if size(z, 1) == meshData.nf
        zCandidates = z(fIdx, :);
        [~, idx] = max(real(conj(oldVel) .* zCandidates), [], 2, 'linear');
        newVel = zCandidates(idx);
        newVel = [real(newVel) imag(newVel)];
        newVel(badpts, :) = NaN;
    else
        zCandidates = z(faces(fIdx, :));
        [~, idx] = max(real(oldVel' * zCandidates), [], 2, 'linear');
        newVel = sum(bary .* zCandidates(idx));
        newVel = [real(newVel) imag(newVel)];
    end
end

end