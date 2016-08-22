function z = swe_invNcdf(p, mu, var)
% Inverse Cumulative Distribution Function (CDF) for univariate Normal
% p: p-value
% mu: mean
% var: variance
% z: equivalent z-score
if nargin < 3
  var = 1;
end
if nargin<2
  mu = 0; 
end
z = -sqrt(2).* (sqrt(var) .* erfcinv(2*p)) + mu;  
end
