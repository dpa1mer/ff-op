function [V,F,b,bc] = flower_closed(n)
%FLOWER Construct a flower with resolution n
%fl is set to false if this is a normal domain, true if it is a domain with
%a single hole in the middle.

if nargin<1
    n = 50;
end

V = {};
F = {};
b = {};
bc = {};

% Parameter domain
theta = linspace(0, 2*pi, n+1)';
thetaM = theta(1:(end-1));
flags = sprintf('-V -q30 -a%0.17f', (2*pi/n)^2);

% Central circle
circVerts = [cos(thetaM), sin(thetaM)];
circEdges = [(1:n)', [(2:n)'; 1]];

% Inner circle with colors
[lV,lF] = triangle((0.5-1e-5)*circVerts,circEdges,[],'Flags',flags);
V{end+1} = lV;
F{end+1} = lF;
b{end+1} = (1:n)';
bc{end+1} = repmat([236,226,240]/255, n, 1);

% Space between inner circle and outer circle
[lV,lF] = triangle([0.5*circVerts; (1-1e-5)*circVerts], ...
    [circEdges; n+circEdges],[0 0],'Flags',flags);
V{end+1} = lV;
F{end+1} = lF;
b{end+1} = (1:(2*n))';
bc{end+1} = [repmat([8,81,156]/255, n, 1); repmat([8,48,10]/255, n, 1)];

% How many petals and petal angle
nPetals = 6;
petalRotAngle = 2*pi / nPetals;
R = [cos(petalRotAngle), -sin(petalRotAngle); ...
    sin(petalRotAngle), cos(petalRotAngle)];

% Up-pointing petal
thetaHalf = theta(1:floor(end/2)+1);
petalVertsOut = ...
    5*[0.2*cos(1.8*(thetaHalf-pi/2)+pi/2), sin(thetaHalf).^2] + [0, 1];
petalVertsIn = ...
    (5-1e-5)*[0.2*cos(1.8*(thetaHalf-pi/2)+pi/2), sin(thetaHalf).^2] + ...
    [0, 1];
petalVerts2Norm = petalVertsIn(:,2) - min(petalVertsIn(:,2));
petalVerts2Norm = petalVerts2Norm / max(petalVerts2Norm);
tinyPetalVertsOut = ...
    2*[0.2*cos(1.8*(thetaHalf-pi/2)+pi/2), sin(thetaHalf).^2] + [0, 1.1];
tinyPetalVertsIn = ...
    (2-1e-5)*[0.2*cos(1.8*(thetaHalf-pi/2)+pi/2), sin(thetaHalf).^2] + ...
    [0, 1.1];
tinyPetalVerts2Norm = tinyPetalVertsIn(:,2) - min(tinyPetalVertsIn(:,2));
tinyPetalVerts2Norm = tinyPetalVerts2Norm / max(tinyPetalVerts2Norm);
petalEdges = [(1:numel(thetaHalf))', [(2:numel(thetaHalf)), 1]'];

% Big outside
squareV = 10*[-1,-1; 1,-1; 1,1; -1,1];
squareE = [1,2; 2,3; 3,4; 4,1];
lineV = [squareV; circVerts];
lineE = [squareE; 4+circEdges];
holes = [0 0];
for i=1:nPetals
    currN = size(lineV,1);
    thesePetals = petalVertsOut * R^(i-1);
    lineV = [lineV; thesePetals];
    lineE = [lineE; currN+petalEdges];
    holes = [holes; mean(thesePetals)];
end
[lV,lF] = triangle(lineV, lineE, holes,'Flags',flags);
V{end+1} = lV;
F{end+1} = lF;
b{end+1} = unique(outline(lF));
normb = normrow(lV(b{end},:));
mnb = min(normb);
mxb = max(normrow(petalVertsOut));
lerp = min((normb-mnb)/(mxb-mnb), 1);
bc{end+1} = lerp .* ([247,252,253]/255) + (1-lerp) .* ([161,217,155]/255);

% Do every petal
for i=1:nPetals
    %Complete inside of petal
    [lV,lF] = triangle(tinyPetalVertsIn * R^(i-1), petalEdges, ...
        [],'Flags',flags);
    V{end+1} = lV;
    F{end+1} = lF;
    b{end+1} = (1:size(tinyPetalVertsIn,1))';
    oColor = tinyPetalVerts2Norm.*[240,240,240]/255 + ...
        (1-tinyPetalVerts2Norm).*[247,252,185]/255;
    bc{end+1} = oColor;
    
    %Space between inside and outside
    [lV,lF] = triangle([petalVertsIn;tinyPetalVertsOut] * R^(i-1), ...
        [petalEdges; size(petalVertsIn,1)+petalEdges], ...
        mean(tinyPetalVertsOut*R^(i-1)),'Flags',flags);
    V{end+1} = lV;
    F{end+1} = lF;
    b{end+1} = [(1:size(petalVertsIn,1))'; ...
        size(petalVertsIn,1)+(1:size(tinyPetalVertsOut,1))'];
    oColor1 = (1-petalVerts2Norm).*[8,48,107]/255 + ...
        petalVerts2Norm.*[43,140,190]/255;
    oColor2 = (1-tinyPetalVerts2Norm).*[4,90,141]/255 + ...
        tinyPetalVerts2Norm.*[247,252,185]/255;
        %tinyPetalVerts2Norm.*[78,179,211]/255;
    bc{end+1} = [oColor1; oColor2];
end

end

