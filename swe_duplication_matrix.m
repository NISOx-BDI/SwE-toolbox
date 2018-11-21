function d = swe_duplication_matrix(n)
% Compute the duplication matrix of size n^2 X n(n+1)/2
% =========================================================================
% FORMAT: d = swe_duplication_matrix(n)
% -------------------------------------------------------------------------
% Inputs:
%  - n: size parameter
% -------------------------------------------------------------------------
% Outputs:
%  - d: duplication matrix
% =========================================================================
% By Bryan Guillaume
% Version Info:  $Format:%ci$ $Format:%h$

    d = zeros(n^2, n * (n + 1) / 2);
    it = 0;
    for j = 1 : n
        d ((j - 1) * n + j, it + j) = 1;
        for i = (j + 1) : n
            d ((j - 1) * n + i, it + i) = 1;
            d ((i - 1) * n + j, it + i) = 1;
        end
        it = it + n - j;
    end
end