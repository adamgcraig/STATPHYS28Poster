function power_spectra = get_power_spectra_all_areas(time_series, sampling_frequency)
%UNTITLED19 Summary of this function goes here
%   Detailed explanation goes here

if ~exist('sampling_frequency', 'var')
    sampling_frequency = 1;
end
num_areas = size(time_series,1);
power_spectra_cell = cell(num_areas,1);
for area_index = 1:num_areas
    power_spectra_cell{area_index} = get_power_spectrum( time_series(area_index,:), sampling_frequency );
end
power_spectra = vertcat(power_spectra_cell{:});

end