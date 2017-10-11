function d = empirical_dist(norm_type,x,y)

n = size(x,2);

% wasserstein integer conversion
Sc = 1e+6;

switch norm_type
    case {'wass1' 'wass2'}
        p = str2num(norm_type(5));
        C0 = distmat(x,y).^p;
        C = int32( Sc*C0 );
        [rho,varrho,u,v] = hungarianLSAP(C);
        d = sum( C0( (1:n) + (double(rho)'-1)*n ) ) / n;
        d = d^(1/p);
    case 'energy'
        k = @(a)-abs(a);
        mu = ones(n,1)/n;
        d = sqrt( rkhs_norm_invariant(k,x,mu,y,mu) );
    case 'gaussian'
        sigma = .3;
        k = @(a)exp( -(a).^2 / (2*sigma^2) );
        mu = ones(n,1)/n;
        d = sqrt( rkhs_norm_invariant(k,x,mu,y,mu) );
    otherwise
        error('Unknown');
end

end