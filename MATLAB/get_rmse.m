function rmse = get_rmse(matrix_1, matrix_2, dim)
%GET_RMSE Compute the root mean squared error between two matrices.
%   matrix_1 and matrix_2 must be the same size in all non-singleton dims.
%   dim gets passed along to mean.
%   If not specified, we use 'all' as the default instead of 1.

if ~exist('dimension', 'var')
    dim = 'all';
end
rmse = sqrt(  mean( (matrix_2 - matrix_1).^2, dim )  );

end