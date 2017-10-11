function [PD, PDinf] = powerDiagramWrapper(E, wts)
% function [PD, PDinf] = powerDiagramWrapper(E, wts)
%
% E: set of points
% wts: weights of the points in E
%
% The output cell PD contains the pieces of the power diagram indexed by
% dimension. PD{1} contains fully-dimensional regions of the power diagram.
% PDinf contains points on infinite edges of the power diagram.
%
% If the initial points are in R^2, the power diagram is drawn.

disp('Lifting points');
LE = liftPD(E, wts);
disp('Computing convex hull');
C = convhulln(LE);
disp('Extracting lower hull');
ind = normalsPD(LE, C);
T = C(ind, :);

nT = size(T,2);
disp('Finding pieces');
[P, total] = piecesPD(T);
disp('Finding power centers');
[PC, powers] = powercentersPD(T, E, wts);
disp('Finding free boundary');
FF = freeBouPD(T, P{nT-1});
disp('Finding power diagram');
PD = pwrDiagramPD(T, PC);

center = mean(E,1);
PDinf = zeros(size(FF));

% find distance from center to farthest powercenter
length = max(sqrt(sum(bsxfun(@minus, PC, center).^2,2)));

disp('Finding points on infinite edges of the power diagram');
for i=1:size(FF,1)
    facet = E(FF(i,:),:);
    ea = edgeAttPD(T, FF(i,:));
    pc = PC(ea{1},:);
    ct = mean(facet,1);
    
    % find vector normal to the facet
    v = null(bsxfun(@minus, facet(1,:), facet(2:end,:)))';
    
    % reorient v to point outward
    if dot(center - ct, v) > 0
        v = -1*v;
    end
    
    % scale v to ensure newpt is sufficiently far away
    v = length*v;
    
    % find point on infinite edge of power diagram
    newpt = pc + v;
    p = piecesPD(FF(i,:));
    for j=1:size(p,1)
        for k=1:size(p{j},1)
            ind = find(ismember(P{j},p{j}(k,:), 'rows'));
            PD{j}{ind} = [PD{j}{ind}; newpt];
        end
    end
    
    % keep track of generated point
    PDinf(i,:) = newpt;
end

% plot 2D power diagram
if (nT-1==2)
    figure
    axis equal
    hold on
    plot(PC(:,1), PC(:,2), 'r.');
    powerdiagram(E(:,1), E(:,2), T, 'r', wts);
    plot(E(:,1), E(:,2), 'b.');
end