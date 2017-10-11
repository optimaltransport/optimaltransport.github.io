function y = inverseFunc(f, maxVal)

% y = inverseFunc(f, maxVal)
% given a monotonous function f: X->[0, maxVal],
% return its inverse f^{-1}
% faster than interp1(f,x)
%author: Nicolas Bonneel

N = length(f);
y = zeros(N, 1);
cursor = 1; 
for i=1:N
    desiredF = i*maxVal/N;
    
    while (cursor<N && f(cursor)<=desiredF)
     cursor = cursor+1;
    end;    
    
    if (cursor<2)
        y(i) = cursor/maxVal;
    else
        alpha = min(1, max(0,(desiredF - f(cursor-1))/(f(cursor)-f(cursor-1))));
        y(i) = ((cursor-1)*(1-alpha) + alpha*cursor )/maxVal;
    end;
end;
