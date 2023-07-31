function [] = save_time_series(ts,file_name)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

file_id = fopen(file_name, 'w');
fwrite(file_id, ts, 'float64');
fclose(file_id);

end