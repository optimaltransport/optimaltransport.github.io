%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is to clip 2 polygons with at least one of them is convex.
% The algorithm uses Sutherland-Hodgman algorithm.
% Conditions:
%   * The clipping polygon must be convex and the subject polygon can be non-convex.
%   * All coordinates must be ordered in CW direction
%   * The last point is not the first point
% Input:
%   * xc: array that specifies the x-coordinate of the clipping polygon (px1)
%   * yc: array that specifies the y-coordinate of the clipping polygon (px1)
%   * xs: array that specifies the x-coordinate of the subject polygon (px1)
%   * ys: array that specifies the y-coordinate of the subject polygon (px1)
% Output:
%   * xp: array that specifies the x-coordinate of the resulted polygon (px1)
%   * yp: array that specifies the y-coordinate of the resulted polygon (px1)
% Extended:
%   * can include jacobian vector, so xc, yc, xs, ys, and xp, yp will have dimension of (px3)
%   * the rows are: x : [xA, dxA/dxi, dxA/dyi; xB, dxB/dxi, dxB/dyi; ...]
%   * the rows are: y : [yA, dyA/dxi, dyA/dyi; yB, dyB/dxi, dyB/dyi, ...]

function [xp, yp] = clip_polygons(xc, yc, xs, ys)
    if (size(xc,2) == 1)
        ixpos = 1;
        iypos = 2;
    else
        ixpos = 1;
        iypos = 4;
    end;
    
    % make the coordinates (x,y) for each polygon
    % if with jacobian, it becomes: [xA, dxA/dxi, dxA/dyi, yA, dyA/dxi, dyA/dyi; ...]
    clippingPolygon = [xc, yc];
    subjectPolygon = [xs, ys];
    
    outputList = subjectPolygon;
    for (i = [1:size(clippingPolygon,1)])
        % break if there are no point left
        if (length(outputList) == 0) break; end;
        
        % get the edge of the clipping polygon
        if (i == size(clippingPolygon,1)) ip1 = 1;
        else ip1 = i+1; end;
        clipEdge = [clippingPolygon(i,:); clippingPolygon(ip1,:) - clippingPolygon(i,:)];
        
        % get the vector pointing inside
        insideVector = [clipEdge(2,iypos), -clipEdge(2,ixpos)];
        
        % copy the output list and clear it
        inputList = outputList;
        outputList = [];
        
        S = inputList(end,:);
        for (iE = [1:size(inputList,1)])
            E = inputList(iE,:);
            SEedge = [S; S-E];
            
            % check if E is inside the clipEdge
            if (isInside(E, clipEdge, insideVector, ixpos, iypos))
            % if (dot(E([ixpos,iypos]) - clipEdge(1,[ixpos,iypos]), insideVector >= 0)
                
                % check if S is not inside the clipEdge
                if (~isInside(S, clipEdge, insideVector, ixpos, iypos))
                % if (dot(S([ixpos,iypos]) - clipEdge(1,[ixpos,iypos]), insideVector) < 0)
                    
                    % add the intersection from S to E with the clipEdge
                    outputList(end+1,:) = getIntersection(SEedge, clipEdge, ixpos, iypos);
                    
                    % A = [SEedge(2,iypos), -SEedge(2,ixpos); clipEdge(2,iypos), -clipEdge(2,ixpos)];
                    % b = [det(SEedge); det(clipEdge)];
                    % intersection = A \ b;
                    % outputList(end+1,:) = intersection';
                end
                outputList(end+1,:) = E;
            
            % check if S is inside the clipEdge
            elseif (isInside(S, clipEdge, insideVector, ixpos, iypos))
            % elseif (dot(S([ixpos,iypos]) - clipEdge(1,[ixpos,iypos]), insideVector) >= 0)
                
                % add the intersection from S to E with the clipEdge
                outputList(end+1,:) = getIntersection(SEedge, clipEdge, ixpos, iypos);
                
                % A = [SEedge(2,iypos), -SEedge(2,ixpos); clipEdge(2,iypos), -clipEdge(2,ixpos)];
                % b = [det(SEedge); det(clipEdge)];
                % intersection = A \ b;
                % outputList(end+1,:) = intersection';
            end
            
            S = E;
        end
    end
    
    if (length(outputList) == 0)
        xp = [];
        yp = [];
    else
        xp = outputList(:,ixpos:iypos-1);
        yp = outputList(:,iypos:end);
    end
end

function ret = isInside(E, clipEdge, insideVector, ixpos, iypos)
    ret = (dot(E([ixpos,iypos]) - clipEdge(1,[ixpos,iypos]), insideVector) >= 0);
end

function intersection = getIntersection(SEedge, clipEdge, ixpos, iypos)
    A = [SEedge(2,iypos), -SEedge(2,ixpos); clipEdge(2,iypos), -clipEdge(2,ixpos)];
    b = [det2(SEedge(:,[ixpos,iypos])); det2(clipEdge(:,[ixpos,iypos]))];
    xy = solve(A, b);
    
    if (iypos == 4)
        Ax = [SEedge(2,iypos+1), -SEedge(2,ixpos+1); clipEdge(2,iypos+1), -clipEdge(2,ixpos+1)];
        bx = [det2(SEedge([1,7;4,10])) + det2(SEedge([3,9;2,8])); det2(clipEdge([1,7;4,10])) + det2(clipEdge([3,9;2,8]))];
        dxydx = solve(A, (bx - Ax*xy));
        
        Ay = [SEedge(2,iypos+2), -SEedge(2,ixpos+2); clipEdge(2,iypos+2), -clipEdge(2,ixpos+2)];
        by = [det2(SEedge([1,7;6,12])) + det2(SEedge([5,11;2,8])); det2(clipEdge([1,7;6,12])) + det2(clipEdge([5,11;2,8]))];
        dxydy = solve(A, (by - Ay*xy));
        intersection = [xy(1) dxydx(1) dxydy(1) xy(2) dxydx(2) dxydy(2)];
    else
        intersection = xy';
    end
end

function r = det2(A)
    r = A(1)*A(4) - A(2)*A(3);
end

function r = solve(A,b)
    r = [ A(4)*b(1)-A(3)*b(2) ; A(1)*b(2)-A(2)*b(1) ] / ( A(1)*A(4) - A(2)*A(3) );
end
