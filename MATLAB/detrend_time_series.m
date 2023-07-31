function ts_detrended = detrend_time_series(ts)
%DETREND_TIME_SERIES Use linear regression to detrend the time series.
%   Do linear regression for each area (row) separately.

ts_detrended = NaN( size(ts) );
num_time_points = size(ts,2);
x = [ ones(num_time_points,1) (1:num_time_points)' ];
for area = 1:size(ts,1)
    y = ts(area,:)';
    p = x\y;
    y_pred = x*p;
    % plot( x(:,2), y, '-r', x(:,2), y_pred, '--g' )
    ts_detrended(area,:) = ( y - y_pred )';
end

end