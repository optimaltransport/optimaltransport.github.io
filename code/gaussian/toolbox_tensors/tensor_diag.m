function D = tensor_diag(a,b)

% diagonal tensor creation

resh = @(x)reshape(x, [1 1 size(x)]);
D = zeros(2,2,size(a,1),size(a,2));
a = resh(a); 
b = resh(b); 
D(1,1,:,:) = a; D(2,2,:,:) = b;

end