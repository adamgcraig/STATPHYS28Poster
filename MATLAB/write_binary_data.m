function [] = write_binary_data(data, file_name)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

file_id = fopen(file_name, 'w');
fwrite(file_id, data, 'float64');
fclose(file_id);

end