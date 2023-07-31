function power_spectrum = get_power_spectrum(time_series, sampling_frequency)
%GET_POWER_SPECTRUM Compute the power spectrum of a time series.
%   We assume the time series
%   is 1-dimensional,
%   is real-valued,
%   and has an even number of samples.
%   See the first example in
%   https://ww2.mathworks.cn/help/signal/ug/power-spectral-density-estimates-using-fft.html

ts_length = numel(time_series);
if ~exist('sampling_frequency','var')
    sampling_frequency = ts_length;
end
ts_dft = fft(time_series);
ts_dft = ts_dft(1:ts_length/2+1);
power_spectrum = ( 1/(ts_length*sampling_frequency) ) * abs(ts_dft).^2;
power_spectrum(2:end-1) = 2 * power_spectrum(2:end-1);

end