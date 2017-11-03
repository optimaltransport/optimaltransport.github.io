function p = mat2perm(S)

p = [];
for i=1:size(S,1)
    [~,p(i)] = max(S(i,:));
end

end