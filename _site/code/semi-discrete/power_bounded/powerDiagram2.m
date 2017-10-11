% This function obtains the power diagram specified by sites E with weights wts.
% It is optimised to work fast on large number of sites (e.g. 10000 sites or more).
% Only works for 2 dimensions.
% Input:
%   * E: a matrix that specifies the sites coordinates (Npts x 2)
%   * wts: a column vector that specifies the sites' weights (Npts x 1)
% Output:
%   * V: list of points' coordinates that makes the power diagram vertices.
%   * CE: cells that contains index of coordinate in V that makes the power diagram of a specified site.
% Thanks to: Arlind Nocaj and Ulrik Brandes (http://onlinelibrary.wiley.com/doi/10.1111/j.1467-8659.2012.03078.x/pdf)
%            and Frederick McCollum (http://uk.mathworks.com/matlabcentral/fileexchange/44385-power-diagrams)
% Made by: Muhammad Kasim (muhammad.kasim@wolfson.ox.ac.uk)

function [V, CE] = powerDiagram2(E, wts)
    %%%%%%%%%%%%%%%%%%%% lift the sites and extend to 3 dimensions %%%%%%%%%%%%%%%%%%%%
    E = [E, sum(E.^2,2)-wts];
    
    %%%%%%%%%%%%%%%%%%%% get the convex hull index %%%%%%%%%%%%%%%%%%%%
    C = convhulln(E);
    
    %%%%%%%%%%%%%%%%%%%% get the lower hull %%%%%%%%%%%%%%%%%%%%
    % centre of the convex hull
    centre = mean(E, 1);
    
    % get the normals
    vec1 = zeros([size(C,1), size(E,2)]);
    vec2 = zeros([size(C,1), size(E,2)]);
    for (i = [1:size(E,2)])
        EiC = E(C+(i-1)*size(E,1));
        vec1(:,i) = EiC(:,2) - EiC(:,1);
        vec2(:,i) = EiC(:,3) - EiC(:,1);
    end
    normals = cross(vec1, vec2, 2);
    vec1 = []; vec2 = []; % remove the memory
    
    % get the middle point of each facet 
    midPoints = zeros([size(C,1), size(E,2)]);
    for (i = [1:size(E,2)])
        EiC = E(C+(i-1)*size(E,1));
        midPoints(:,i) = mean(EiC,2);
    end
    EiC = []; % remove the memory
    
    % check the projections of normals to midPoints-centre vector
    dot = sum(bsxfun(@minus, centre, midPoints) .* normals, 2);
    outward = (dot < 0);
    
    % make sure all normals are pointing inwards
    normals(outward,:) = -normals(outward,:);
    
    % get the lower hull & upper hull
    lowerIdx = (normals(:,end) > 0);
    upperIdx = (normals(:,end) <= 0);
    CUp = C(upperIdx,:);
    C = C(lowerIdx,:);
    normals = normals(lowerIdx,:);
    midPoints = midPoints(lowerIdx,:);
    
    %%%%%%%%%%%%%%%%%%%% invert the facet from dual space to the real space %%%%%%%%%%%%%%%%%%%%
    normalsZ = bsxfun(@rdivide, normals, -normals(:,end)); % normalise the normals to have normals(z) = -1
    a = normalsZ(:,1);
    b = normalsZ(:,2);
    % c = midPoints(:,3) - a.*midPoints(:,1) - b.*midPoints(:,2);
    V = [a/2, b/2]; % this is the vertices for the power diagrams
    V = [Inf, Inf; V];
    
    %%%%%%%%%%%%%%%%%%%% assign each point to the vertices %%%%%%%%%%%%%%%%%%%%
    % CE = arrayfun(@(x) mod(find(C == x)'-1, size(C,1)) + 1, [1:size(E,1)], 'UniformOutput', 0);
    CE = cell(size(E,1),1);
    for (col = [1:size(C,2)])
        for (row = [1:size(C,1)])
            i = C(row,col);
            CE{i} = [CE{i}, row+1]; % + 1 because there is Inf at the first row
        end
    end
    
    % select which one is on border
    onBorders = zeros([size(E,1),1]);
    for (col = [1:size(CUp,2)])
        for (row = [1:size(CUp,1)])
            i = CUp(row,col);
            if (onBorders(i)) continue; end;
            onBorders(i) = 1;
            CE{i} = [CE{i}, 1];
        end
    end
end
