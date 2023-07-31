% PLOT_SINGLE_TS_LINEAR_MODEL_WEIGHT_VARIABILITY
% Adam Craig, 2023-06-09
% Visualize the variability of weights in single-ts linear models.

hcp_data_header

print_status_update_if_time('starting code for comparison of single-ts linear model weights...')

use_y_intercept = true;
use_y_intercept_string = sprintf('yint_%u', use_y_intercept);

use_sc = false;
use_sc_string = sprintf('use_sc_%u', use_sc);

% nonlinearity = @(v) (2/pi)*atan(v);
% nonlinearity_inverse = @(v) tan( (pi/2)*v );
% nonlinearity_name = 'nl_atan';
% rescale_name = 'min_max_norm';
% rescale_fun = @(ts) rescale_ts( ts, -0.999, 0.999 );
% rand_min = -0.999;
% rand_max = 0.999;

% nonlinearity = @(v) tanh(v);
% nonlinearity_inverse = @(v) atanh( v );
% nonlinearity_name = 'nl_tanh';
% rescale_name = 'min_max_norm';
% % Make the range a little smaller than [-1, +1]
% % so that atanh does not give any infinite values.
% rescale_fun = @(ts) rescale_ts( ts, -0.999, 0.999 );
% rand_min = -0.999;
% rand_max = 0.999;

% nonlinearity = @(v) 1./( 1 + exp(-v) );
% nonlinearity_inverse = @(v) -log( 1./v - 1 );
% nonlinearity_name = 'nl_logistic';
% rescale_name = 'min_max_norm_positive';
% rescale_fun = @(ts) rescale_ts( ts, 0.001, 0.999 );
% rand_min = 0.001;
% rand_max = 0.999;

% leak_size = 0.001;
% nonlinearity = @(v) ( (v >= 0) + (v < 0)*leak_size ).*v;
% nonlinearity_inverse = @(v) ( (v >= 0) + (v < 0)/leak_size ).*v;
% nonlinearity_name = 'nl_leaky_relu';
% rescale_name = 'min_max_norm_positive';
% rescale_fun = @(ts) rescale_ts( ts, 0.001, 0.999 );
% rand_min = 0.001;
% rand_max = 0.999;

% nonlinearity = @(v) v;
% nonlinearity_inverse = @(v) v;
% nonlinearity_name = 'raw';
% rescale_name = 'raw';
% rescale_fun = @(ts) ts;
% rand_min = -10.0;
% rand_max = 10.0;

% nonlinearity = @(v) v;
% nonlinearity_inverse = @(v) v;
% nonlinearity_name = 'zero_mean';
% rescale_name = 'zero_mean';
% rescale_fun = @(ts) ts - mean(ts,'all');
% rand_min = -10.0;
% rand_max = 10.0;

% nonlinearity = @(v) v;
% nonlinearity_inverse = @(v) v;
% nonlinearity_name = 'detrend';
% rescale_name = 'detrend';
% rescale_fun = @(ts) detrend_time_series(ts);
% rand_min = -10.0;
% rand_max = 10.0;

% nonlinearity = @(v) v;
% nonlinearity_inverse = @(v) v;
% nonlinearity_name = 'std_mean_norm';
% rescale_name = 'std_mean_norm';
% rescale_fun = @(ts) std_mean_normalize_ts(ts);
% rand_min = -2.0;
% rand_max = 2.0;

% nonlinearity = @(v) v;
% nonlinearity_inverse = @(v) v;
% nonlinearity_name = 'min_max_norm_positive';
% rescale_name = 'min_max_norm_positive';
% rescale_fun = @(ts) rescale_ts( ts, 0.0, 1.0 );
% rand_min = 0.0;
% rand_max = 1.0;

nonlinearity = @(v) v;
nonlinearity_inverse = @(v) v;
nonlinearity_name = 'min_max_norm';
rescale_name = 'min_max_norm';
rescale_fun = @(ts) rescale_ts( ts, -1.0, 1.0 );
rand_min = -1.0;
rand_max = 1.0;

settings_string = [ use_y_intercept_string '_' use_sc_string '_' nonlinearity_name ];

group = 'training';
subject_ids = training_subject_ids;
% group = 'validation';
% subject_ids = validation_subject_ids;
% group = 'testing';
% subject_ids = testing_subject_ids;

num_subjects = numel(subject_ids);

W = NaN(num_brain_areas, num_brain_areas+1, num_subjects * time_series_per_subject);
total_ts_index = 1;
for subject_index = 1:num_subjects
    subject_id = subject_ids(subject_index);
    for time_series_index = 1:time_series_per_subject
        ts_suffix = time_series_strings{time_series_index};
        model_file_name = [single_ts_linear_model_dir sprintf('model_%s_%u_%s.bin', settings_string, subject_id, ts_suffix)];
        W(:,:,total_ts_index) = load_data_from_binary(model_file_name, num_brain_areas, num_brain_areas+1);
        print_status_update_if_time( sprintf('subject %u of %u, time series %u of %u', ...
            subject_index, num_subjects, time_series_index, time_series_per_subject) )
        total_ts_index = total_ts_index + 1;
    end
end

fig_mean = figure;
W_mean = mean(W,3);
heatmap(W_mean,'Title','mean weight','GridVisible','off')

fig_std = figure;
W_std = std(W,0,3);
heatmap(W_std,'Title','std. dev. of weight','GridVisible','off')

fig_min = figure;
W_min = min(W,[],3);
heatmap(W_min,'Title','minimum weight','GridVisible','off')

fig_max = figure;
W_max = max(W,[],3);
heatmap(W_max,'Title','maximum weight','GridVisible','off')

fig_range = figure;
W_range = W_max - W_min;
heatmap(W_range,'Title','range of weight','GridVisible','off')

fig_b_hist = figure;
histogram(  reshape( W(:,1,:), [], 1, 1 ), 'Normalization', 'probability'  )
xlabel('linear model bias term for area')
ylabel('probability')

fig_W_hist = figure;
histogram(  reshape( W(:,2:end,:), [], 1, 1 ), 'Normalization', 'probability'  )
xlabel('linear model weight between pairs of areas')
ylabel('probability')

fig_b_range_hist = figure;
histogram(  reshape( W_range(:,1), [], 1 ), 'Normalization', 'probability'  )
xlabel('range of linear model bias term for area')
ylabel('probability')

fig_W_range_hist = figure;
histogram(  reshape( W_range(:,2:end), [], 1 ), 'Normalization', 'probability'  )
xlabel('range of linear model weight between pairs of areas')
ylabel('probability')

fig_b_std_hist = figure;
histogram(  reshape( W_std(:,1), [], 1 ), 'Normalization', 'probability'  )
xlabel('std. dev. of linear model bias term for area')
ylabel('probability')

fig_W_std_hist = figure;
histogram(  reshape( W_std(:,2:end), [], 1 ), 'Normalization', 'probability'  )
xlabel('std. dev. of linear model weight between pair of areas')
ylabel('probability')

