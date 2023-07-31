function data_matrix = load_data_from_binary(binary_data_file, num_rows, num_cols)
%LOAD_DATA_FROM_BINARY Load data from a binary file.
%   We need to specify the dimensions,
%   because we store the data as an array of 64-bit floats in a binary.

file_id = fopen(binary_data_file,'r');
binary_data = fread(file_id, num_rows*num_cols, 'float64');
fclose(file_id);
data_matrix = reshape(binary_data, num_rows, num_cols);

end