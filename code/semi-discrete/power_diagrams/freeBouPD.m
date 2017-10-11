function FB = freeBouPD(T, P)
% function FB = freeBouPD(T, P)
%
% T: triangulation
% P: pieces of the triangulation
%
% The output FB contains the pieces of P on the boundary of T.

[~, nT] = size(T);
[mP, ~] = size(P);

ii=1;
for i=1:mP
    if size(find(sum(ismember(T, P(i,:)),2)==nT-1),1)==1
        FB(ii,:) = P(i,:);
        ii = ii+1;
    end
end