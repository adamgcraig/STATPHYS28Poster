function before_and_after_table = load_time_series_as_before_and_after_table(ts_file_name)
%LOAD_TIME_SERIES_AS_BEFORE_AND_AFTER_TABLE Load ts as before-after table.
%   ts_file_name: time series loadable with load_time_series()
%   before_and_after_table: 2 columns, before and after.

data_ts = load_time_series(ts_file_name)';
before = data_ts(1:end-1,:);
after = data_ts(2:end,:);
before_and_after_table = table(before, after);


end