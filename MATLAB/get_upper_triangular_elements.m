function upper_triangular_elements = get_upper_triangular_elements(two_d_matrix)
%GET_UPPER_TRIANGULAR_ELEMENTS Get a vector with the upper diagonal elements.
%   By this, we mean the elements above the diagonal,
%   in other words, any two_d_matrix(i,j) where i < j.

[num_rows, num_cols] = size(two_d_matrix);
is_upper_triangular = (1:num_rows)' < 1:num_cols;
upper_triangular_elements = two_d_matrix(is_upper_triangular);

end