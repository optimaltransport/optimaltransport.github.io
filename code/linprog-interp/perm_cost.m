function cost = perm_cost(p,c)

cost = 0;
for q=1:length(p)
    cost = cost + c(q,p(q));
end

end