function C = tensor_mult(A,B)

% tensor_mult - multiply 2 tensors
%
% C = tensor_mult(A,B);
%
%   C = A*B
%
%   Copyright (c) 2016 Gabriel Peyre

[mA,nA,NA1,NA2] = size(A);
[mB,nB,NB1,NB2] = size(B);

if NA1~=NB1 || NA2~=NB2 || nA~=mB
    error('Size problem');
end
C = zeros(mA,nB,NA1,NA2);

for i=1:mA
    for j=1:nB
        for k=1:nA
            C(i,j,:,:) = C(i,j,:,:) + A(i,k,:,:).*B(k,j,:,:);
        end
    end
end

end
