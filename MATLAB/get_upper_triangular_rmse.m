function rmse = get_upper_triangular_rmse( matrix_1, matrix_2 )
%GET_UPPER_TRIANGULAR_RMSE Get RMSE between corresponding upper triangles.
%   rmse is a single scalar value.

rmse = get_rmse( get_upper_triangular_elements(matrix_1), get_upper_triangular_elements(matrix_2), 'all' );

end