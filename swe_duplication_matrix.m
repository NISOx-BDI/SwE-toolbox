function d = swe_duplication_matrix(n)
% Compute the duplication matrix of size n^2 X n(n+1)/2
% By Bryan Guillaume
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
