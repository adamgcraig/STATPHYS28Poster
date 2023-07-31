function x = generate_time_series_with_linear_model(W, x_0, num_steps, use_y_intercept, sc_data, nonlinearity)
%GENERATE_TIME_SERIES_WITH_LINEAR_MODEL Generate time series.
%   x(t+1) = W*x(t)
%   To include a non-0 y-intercept, set use_y_intercept to true.
%   We then require that W is Nx(N+1) instead of NxN.
%   To evaluate the model, we separate out the first column of W:
%   b = W(:,1), W = (:,2:end).
%   Then x(:,t+1) = b + W*x(:,t).
%   By passing in a function handle for nonlinearity,
%   we can use a model x(:,t+1) = f( b + W * x(:,t) )
%   where nonlinearity = f(.).

num_areas = size(x_0, 1);
if ~exist('use_y_intercept', 'var')
    use_y_intercept = false;
end
if ~exist('sc_data', 'var')
    sc_data = zeros(num_areas, num_areas);
end
if ~exist('nonlinearity', 'var')
    nonlinearity = @(v) v;
end

if use_y_intercept
    b = W(:,1);
    W = W(:,2:end);
else
    b = zeros(num_areas, 1);
end

x = NaN(num_areas, num_steps);
x_t = x_0;
I_plus_sc = eye(num_areas, num_areas) + sc_data;
for t = 1:num_steps
    x_t = nonlinearity( b + W*I_plus_sc*x_t );
    x(:,t) = x_t;
end

end