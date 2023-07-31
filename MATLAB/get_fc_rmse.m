function rmse = get_fc_rmse(matrix_1, matrix_2)
%GET_FC_RMSE Compute the root mean squared error between upper triangles.
%   matrix_1 and matrix_2 must be same-size square matrices.

[rows, cols] = size(matrix_1);
ut_indices = (1:rows)' < (1:cols);
rmse = get_rmse( matrix_1(ut_indices), matrix_2(ut_indices), 'all' );

end