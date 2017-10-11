function PD = pwrDiagramPD(T, PC)
% function PD = pwrDiagramPD(T, PC)
%
% T: triangulation
% PC: power centers of triangles in T
%
% The output cell PD contains pieces of the power diagram, indexed by
% dimension. PD{mP} contains pieces of dimension zero (power centers that
% correspond to fully-dimensional pieces of the triangulation T) and PD{1}
% contains fully-dimensional regions of the power diagram (corresponding to
% vertices of the triangulation).

P = piecesPD(T);
mP = size(P,1);
PD = cell(mP,1);

for i=1:mP
    EA = edgeAttPD(T, P{i});
    for j=1:size(EA,1);
        EA{j} = PC(EA{j},:);
    end
    PD{i} = EA;
end