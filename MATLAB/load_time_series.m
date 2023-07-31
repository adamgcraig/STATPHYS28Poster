function time_series = load_time_series(time_series_file, num_brain_areas, num_time_points, offset)
%LOAD_TIME_SERIES Loads an HCP fMRI time series from a binary file.
%   We need to specify the dimensions,
%   because we store the data as an array of 64-bit floats in a binary.

if ~exist('num_brain_areas','var')
    num_brain_areas = 360;
end
if ~exist('num_time_points','var')
    num_time_points = 1200;
end
if ~exist('offset','var')
    offset = 0;
end
time_series = load_data_from_binary(time_series_file, num_brain_areas, offset+num_time_points);
time_series = time_series(:,offset+1:offset+num_time_points);

end