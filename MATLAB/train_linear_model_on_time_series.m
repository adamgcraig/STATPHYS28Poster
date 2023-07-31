function W = train_linear_model_on_time_series(time_series_set, use_y_intercept, sc_data, nonlinearity_inverse)
%TRAIN_LINEAR_MODEL Train a linear model on multiple time series.
%   The linear model is just a matrix W such that
%   prediction for x(:,t+1) = W * x(:,t).
%   To include a non-0 y-intercept, set use_y_intercept to true.
%   We then add a row of 1s to the top of x so that
%   W is Nx(N+1) instead of NxN.
%   To evaluate the model, you can then either prepend a 1 to x
%   or separate out the first column of W:
%   b = W(:,1), W = (:,2:end).
%   Then x(:,t+1) = b + W*x(:,t).
%   By passing in a function handle for nonlinearity_inverse,
%   we can use a model x(:,t+1) = f( W * x(:,t) )
%   where nonlinearity_inverse = f^-1(.) so that we can solve for W
%   by performing a linear regression on
%   f^-1( x(:,t+1) ) vs x(:,t)

num_brain_areas = size(time_series_set{1}, 1);
if ~exist('use_y_intercept', 'var')
    use_y_intercept = false;
end
if ~exist('sc_data', 'var')
    sc_data = zeros(num_brain_areas, num_brain_areas);
end
if ~exist('nonlinearity_inverse', 'var')
    nonlinearity_inverse = @(v) v;
end
x_cell = cellfun( @(ts) ts(:,1:end-1), time_series_set, 'UniformOutput', false);
x = ( eye(num_brain_areas, num_brain_areas) + sc_data ) * horzcat(x_cell{:});
if use_y_intercept
    x = [ ones( 1, size(x,2) ); x ];
end
y_cell = cellfun( @(ts) ts(:,2:end), time_series_set, 'UniformOutput', false);
y = nonlinearity_inverse( horzcat(y_cell{:}) );
W = y/x;

end