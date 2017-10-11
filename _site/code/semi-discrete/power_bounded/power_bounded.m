% POWER_BOUNDED computes the power cells about the points (x,y) inside
% the bounding box (must be a rectangle or a square) crs.  If crs is not supplied, an
% axis-aligned box containing (x,y) is used.
% It is optimised to work fast on large number of sites (e.g. 10000 sites or more)
% Input:
%   * x, y: coordinate of the Voronoi point (numPoints x 1)
%   * wts: weights of each point (numPoints x 1)
%   * crs: vortices of the bounding box in cw order (numVert x 2)
% Output:
%   * V: x,y-coordinate of vertices of the power cells
%   * C: indices of the Voronoi cells from V
% See Matlab's voronoin for more information about the output
% Made by: Aaron Becker, atbecker@uh.edu, and Muhammad Kasim, muhammad.kasim@wolfson.ox.ac.uk

function [V,C] = power_bounded(x,y, wts, crs)
    bnd=[min(x) max(x) min(y) max(y)]; %data bounds
    if nargin < 3
        crs=double([bnd(1) bnd(4);bnd(2) bnd(4);bnd(2) bnd(3);bnd(1) bnd(3);bnd(1) bnd(4)]);
    end

    rgx = max(crs(:,1))-min(crs(:,1));
    rgy = max(crs(:,2))-min(crs(:,2));
    rg = max(rgx,rgy);
    midx = (max(crs(:,1))+min(crs(:,1)))/2;
    midy = (max(crs(:,2))+min(crs(:,2)))/2;

    % add 4 additional edges
    xA = [x; midx + [0;0;-5*rg;+5*rg]];
    yA = [y; midy + [-5*rg;+5*rg;0;0]];
    
    if (all(wts == 0))
        [vi,ci] = voronoin([xA,yA]);
    else
        [vi,ci] = powerDiagram2([xA,yA], [wts;zeros(4,1)]);
    end
    
    % remove the last 4 cells
    C = ci(1:end-4);
    V = vi;
    % use Polybool to crop the cells
    %Polybool for restriction of polygons to domain.
    
    maxX = max(crs(:,1)); minX = min(crs(:,1));
    maxY = max(crs(:,2)); minY = min(crs(:,2));
    for ij=1:length(C)
        % thanks to http://www.mathworks.com/matlabcentral/fileexchange/34428-voronoilimit
        Cij = C{ij};
        if (length(Cij) == 0) continue; end;
        
        % first convert the contour coordinate to clockwise order:
        pts = V(Cij,:);
        K = convhull(pts);
        K = K(end-1:-1:1);
        C{ij} = Cij(K);
        X2 = pts(K,1);
        Y2 = pts(K,2);
        
        % if all points are inside the bounding box, then skip it
        if (all((X2 <= maxX) & (X2 >= minX) & (Y2 <= maxY) & (Y2 >= minY))) continue; end;
        
        [xb, yb] = clip_polygons(crs(:,1),crs(:,2),X2,Y2);
        % xb = xb'; yb = yb';
        ix=nan(1,length(xb));
        for il=1:length(xb)
            if any(V(:,1)==xb(il)) && any(V(:,2)==yb(il))
                ix1=find(V(:,1)==xb(il));
                ix2=find(V(:,2)==yb(il));
                for ib=1:length(ix1)
                    if any(ix1(ib)==ix2)
                        ix(il)=ix1(ib);
                    end
                end
                if isnan(ix(il))==1
                    lv=length(V);
                    V(lv+1,1)=xb(il);
                    V(lv+1,2)=yb(il);
                    ix(il)=lv+1;
                end
            else
                lv=length(V);
                V(lv+1,1)=xb(il);
                V(lv+1,2)=yb(il);
                ix(il)=lv+1;
            end
        end
        C{ij} = ix;
    end
end

