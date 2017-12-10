function CP = copula(gamma)

% copula - Compute copular associated to a joint distribution gamma.
%
% CP = copula(gamma)
%
%   Copyright (c) 2017 Gabriel Peyre

N = size(gamma,1);

mu = sum(gamma,2);
nu = sum(gamma,1)';

% cumulative
cmu = [0;cumsum(mu)];
cnu = [0;cumsum(nu)];
% inverse cumulatives
t = (0:N)'/N;
icmu = interp1(cmu, t, t, 'spline');
icnu = interp1(cnu, t, t, 'spline');
%
Cgamma = cumsum(cumsum(gamma,1),2);
s = linspace(0,1,N); 

icmu = max(min(icmu,1),0);
icnu = max(min(icnu,1),0);

[S2,S1] = meshgrid(s,s);
[T2,T1] = meshgrid(icnu,icmu);
CP = interp2(S2,S1,Cgamma,T2,T1);

end