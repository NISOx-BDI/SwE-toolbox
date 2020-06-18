function z = swe_invNcdf(p, mu, var)
% Inverse Cumulative Distribution Function (CDF) for univariate Normal
% =========================================================================
% FORMAT: z=swe_invNcdf(p, mu, var)
% -------------------------------------------------------------------------
% Inputs:
%   - p: p-value
%   - mu: mean
%   - var: variance
% -------------------------------------------------------------------------
% Outputs:
%   - z: equivalent z-score
% =========================================================================
% Version Info:  $Format:%ci$ $Format:%h$

    if nargin < 3
      var = 1;
    end
    if nargin<2
      mu = 0;
    end
    z = -sqrt(2).* (sqrt(var) .* erfcinv(2*p)) + mu;

end
