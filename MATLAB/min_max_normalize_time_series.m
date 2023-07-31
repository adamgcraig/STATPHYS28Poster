function normed_ts = min_max_normalize_time_series(ts)
%MIN_MAX_NORMALIZE_TIME_SERIES Normalize time series to have range [-1,+1].
%   normed_ts = ( 2*ts - min(ts,[],'all') - max(ts,[],'all') )./range(ts,'all')

ts_min = min(ts,[],'all');
ts_max = max(ts,[],'all');
ts_half_range = (ts_max - ts_min)/2;
ts_midpoint = (ts_max + ts_min)/2;
normed_ts = (ts - ts_midpoint)./ts_half_range;

end