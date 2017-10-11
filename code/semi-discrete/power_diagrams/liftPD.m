function LE = liftPD(E, wts)
% function LE = liftPD(E, wts)
%
% E: set of points
% wts: weights of the points in E
%
% The output array LE contains the points of E lifted one dimension higher.

[N, d] = size(E);
LE = zeros(N,d+1);

for i=1:N
    x = E(i,:);
    LE(i,:) = [x, x*x' - wts(i)];
end