function fc_corr = get_fc_corr( matrix_1, matrix_2 )
%GET_FC_CORR Get correlation between upper triangular elements of matrices.
%   Detailed explanation goes here

[rows, cols] = size(matrix_1);
ut_indices = (1:rows)' < (1:cols);
fc_corr = corr( matrix_1(ut_indices), matrix_2(ut_indices) );

end