function fc = get_functional_connectivity(time_series)
%GET_FUNCTIONAL_CONNECTIVITY Get functional connectivity of time series. 
%   areas are in rows, time points in columns.

fc = corr(time_series');

end