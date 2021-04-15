function [curveVerts,curveEdges,colors] = flower(n)
%FLOWER Construct a flower with resolution n

if nargin<1
    n = 50;
end

% Parameter domain
theta = linspace(0, 2*pi, n+1)';
thetaM = theta(1:(end-1));

% Central circle
circVerts = [cos(thetaM), sin(thetaM)];
circEdges = [(1:n)', [(2:n)'; 1]];
circColIn = repmat([8,48,10]/255, n, 1);
circColOut = repmat([247,252,253]/255, n, 1);
tinyCircColOut = repmat([8,81,156]/255, n, 1);
tinyCircColIn = repmat([236,226,240]/255, n, 1);

curveVerts = {circVerts, (1-1e-5)*circVerts,...
    0.5*circVerts, (0.5-1e-5)*circVerts};
curveEdges = {circEdges, circEdges, circEdges, circEdges};
colors = {circColOut, circColIn, tinyCircColOut, tinyCircColIn};

% Up-pointing petal
thetaHalf = theta(1:floor(end/2)+1);
petalVertsOut = ...
    5*[0.2*cos(1.8*(thetaHalf-pi/2)+pi/2), sin(thetaHalf).^2] + [0, 1];
petalVertsIn = ...
    (5-1e-5)*[0.2*cos(1.8*(thetaHalf-pi/2)+pi/2), sin(thetaHalf).^2] + ...
    [0, 1];
tinyPetalVertsOut = ...
    2*[0.2*cos(1.8*(thetaHalf-pi/2)+pi/2), sin(thetaHalf).^2] + [0, 1.1];
tinyPetalVertsIn = ...
    (2-1e-5)*[0.2*cos(1.8*(thetaHalf-pi/2)+pi/2), sin(thetaHalf).^2] + ...
    [0, 1.1];
petalEdges = [(1:(numel(thetaHalf)-1))', (2:numel(thetaHalf))'];

% How many petals and petal angle
nPetals = 6;
petalRotAngle = 2*pi / nPetals;
R = [cos(petalRotAngle), -sin(petalRotAngle); ...
    sin(petalRotAngle), cos(petalRotAngle)];

for i=1:nPetals
    curveVerts{end+1} = petalVertsOut * R^(i-1);
    curveVerts{end+1} = petalVertsIn * R^(i-1);
    curveEdges{end+1} = petalEdges;
    curveEdges{end+1} = petalEdges;
    colors{end+1} = repmat([247,252,253]/255, size(petalVertsOut,1), 1);
    colors{end+1} = repmat([43,140,190]/255, size(petalVertsOut,1), 1);
    curveVerts{end+1} = tinyPetalVertsOut * R^(i-1);
    curveVerts{end+1} = tinyPetalVertsIn * R^(i-1);
    curveEdges{end+1} = petalEdges;
    curveEdges{end+1} = petalEdges;
    colors{end+1} = repmat([78,179,211]/255, size(petalVertsOut,1), 1);
    colors{end+1} = repmat([247,252,185]/255, size(petalVertsOut,1), 1);
end

end

