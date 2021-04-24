function DistanceExperiment(meshData, source, ellipticity, bcNeumann)

if nargin < 4
    bcNeumann = false;
end

if nargin < 3
    ellipticity = 1e-2;
end

nv = meshData.nv;
verts = meshData.verts;
faces = meshData.faces;
tr = triangulation(faces, verts(:, 1), verts(:, 2));

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

% source = 2567;%4550;%21084;%randi(length(verts), 1);
sqDist = sum((W - W(source, :)).^2 ./ (diag(D).').^2, 2);
dist = sqrt(sqDist);
figure; trisurf(faces, verts(:, 1), verts(:, 2), verts(:, 3), dist, 'EdgeColor', 'none');
view(2); axis image off; shading interp; colormap viridis(20);
hold on; scatter3(verts(source, 1), verts(source, 2), verts(source, 3), 200, 'r.');

%% Compute shortest paths - adapted from Justin's code

np = 400; % number of points moving around
points = verts(meshData.intIdx(randi(length(meshData.intIdx), np, 1)),1:2); % start one path at each random interior vertex

% convenience matrices for divided difference approximation
eps = 1e-6; % epsilon for divided difference gradients
xShift = zeros(np,2); xShift(:,1) = eps;
yShift = zeros(np,2); yShift(:,2) = eps;

% stepsize for gradient descent
stepsize = 0.5 * min(meshData.edgeLengths);%30 * ellipticity;
% number of steps of gradient descent
nsteps = 5000;

paths = nan(np, 2, nsteps);
paths(:, :, 1) = points;

% Gradient descent on squared distance
for i=1:nsteps
    % Divided difference gradient
    dx = (sqdistToSource(points+xShift)-sqdistToSource(points-xShift)) / (2*eps); 
    dy = (sqdistToSource(points+yShift)-sqdistToSource(points-yShift)) / (2*eps);
    gradient = [dx dy];
    
    % Normalize gradient
    gradient = gradient ./ vecnorm(gradient, 2, 2);
    
    % Gradient descent
    points = points - stepsize*gradient; 
    
    % Add points to paths
    paths(:, :, i) = points;
end

plot(squeeze(paths(:, 1, :))', squeeze(paths(:, 2, :))', 'r', 'LineWidth', 1);

% takes an np x 2 matrix and gives back distance of each of those to "source"
function sqdPts = sqdistToSource(pts)
    % Get barycentric coords
    badpts = any(isnan(pts), 2);
    fIdx = ones(size(pts, 1), 1);
    bary = zeros(size(pts, 1), 3);
    [fIdx(~badpts), bary(~badpts, :)] = pointLocation(tr, pts(~badpts, :));
    badpts = badpts | isnan(fIdx);
    fIdx(badpts) = 1; % placeholder index
    sqdPts = sum(bary .* sqDist(faces(fIdx, :)), 2);
    sqdPts(badpts) = NaN;
end

end