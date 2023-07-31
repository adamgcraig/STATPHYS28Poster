% MAKE_AND_SAVE_SINGLE_TS_LINEAR_MODELS
% Adam Craig, 2023-06-06
% Make and save simple linear models for all the 

hcp_data_header

print_status_update_if_time('starting code for creation of single-ts linear models...')

use_y_intercept = true;
use_y_intercept_string = sprintf('yint_%u', use_y_intercept);

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

nonlinearity = @(v) v;
nonlinearity_inverse = @(v) v;
nonlinearity_name = 'min_max_norm';
rescale_name = 'min_max_norm';
rescale_fun = @(ts) rescale_ts( ts, -1.0, 1.0 );
rand_min = -1.0;
rand_max = 1.0;

% nonlinearity = @(v) v;
% nonlinearity_inverse = @(v) v;
% nonlinearity_name = 'min_max_norm_positive';
% rescale_name = 'min_max_norm_positive';
% rescale_fun = @(ts) rescale_ts( ts, 0.0, 1.0 );
% rand_min = 0.0;
% rand_max = 1.0;

use_sc = false;
use_sc_string = sprintf('use_sc_%u', use_sc);

settings_string = [ use_y_intercept_string '_' use_sc_string '_' nonlinearity_name ];

% group = 'training';
% subject_ids = training_subject_ids;
group = 'validation';
subject_ids = validation_subject_ids;
% group = 'testing';
% subject_ids = testing_subject_ids;

num_subjects = numel(subject_ids);

ts_nans = NaN(num_subjects, time_series_per_subject);
ts_rmse = ts_nans;
ps_rmse = ts_nans;
fc_rmse = ts_nans;
has_sc = false(num_subjects,1);
for subject_index = 1:num_subjects
    subject_id = subject_ids(subject_index);
    sc_file = get_structural_connectivity_file(subject_id);
    if exist(sc_file,'file')
        current_subject_has_sc = true;
        has_sc(subject_index) = true;
        data_sc = load_structural_connectivity(sc_file);
    else
        data_sc = zeros(num_brain_areas, num_brain_areas);
    end
    if use_sc
        sc_input = data_sc;
    else
        sc_input = zeros( size(data_sc) );
    end
    for time_series_index = 1:time_series_per_subject
        ts_suffix = time_series_strings{time_series_index};
        time_series_file = get_time_series_file(subject_id, ts_suffix);
        data_ts = rescale_fun( load_time_series(time_series_file) );
        data_ps = get_power_spectra_all_areas(data_ts, sampling_frequency);
        data_fc = get_functional_connectivity(data_ts);
        W = train_linear_model_on_time_series({data_ts}, use_y_intercept, sc_input, nonlinearity_inverse);
        model_file_name = [single_ts_linear_model_dir sprintf('model_%s_%u_%s.bin', settings_string, subject_id, ts_suffix)];
        write_binary_data(W, model_file_name)
        sim_ts = generate_time_series_with_linear_model( W, data_ts(:,1), num_time_points, use_y_intercept, sc_input, nonlinearity );
        sim_ps = get_power_spectra_all_areas(sim_ts, sampling_frequency);
        sim_fc = get_functional_connectivity(sim_ts);
        ts_rmse(subject_index, time_series_index) = get_rmse( sim_ts(:,2:end), data_ts(:,2:end) );
        ps_rmse(subject_index, time_series_index) = get_rmse(sim_ps, data_ps);
        fc_rmse(subject_index, time_series_index) = get_upper_triangular_rmse(sim_fc, data_fc);
        print_status_update_if_time( sprintf('subject %u of %u, time series %u of %u', ...
            subject_index, num_subjects, time_series_index, time_series_per_subject) )
    end
end
last_subject_id_str = sprintf('subject_%u', subject_id);

% Combine the stored results into a single table.
subject_id = reshape( repmat(subject_ids, 1, time_series_per_subject), [], 1 );
time_series = reshape( repmat(time_series_strings, num_subjects, 1), [], 1 );
has_sc = reshape( repmat(has_sc, 1, time_series_per_subject), [], 1 );
ts_rmse = ts_rmse(:);
ps_rmse = ps_rmse(:);
fc_rmse = fc_rmse(:);
results_table = table(subject_id, time_series, has_sc, ts_rmse, ps_rmse, fc_rmse);
writetable(results_table, ['single_ts_linear_model_' settings_string '_rmses.csv'])

% Compare the FC and PS RMSE distributions to
% the distribution of RMSEs between FCs/PSs
% from different time series of the same subject.
% For comparison, compare the distribution of intra-subject RMSEs to
% the distribution of inter-subject RMSEs.

same_subject_ps_rmse = readmatrix([rmse_matrices_dir group '_fmri_power_spectrum_rmses_same_subject_' rescale_name '.csv']);
different_subject_ps_rmse = readmatrix([rmse_matrices_dir group '_fmri_power_spectrum_rmses_inter_subject_' rescale_name '.csv']);
same_subject_fc_rmse = readmatrix([rmse_matrices_dir group '_fmri_functional_connectivity_rmses_same_subject_' rescale_name '.csv']);
different_subject_fc_rmse = readmatrix([rmse_matrices_dir group '_fmri_functional_connectivity_rmses_inter_subject_' rescale_name '.csv']);

% Do statistical tests to
% see which RMSE distributions are different from each other.
% From https://ww2.mathworks.cn/help/stats/kstest2.html#namevaluepairarguments
% 'Tail', 'smaller'
% "Test the alternative hypothesis that the empirical cdf of x1 is smaller than the empirical cdf of x2."
% "If the data values in x1 tend to be larger than those in x2,
% the empirical distribution function of x1 tends to be smaller than that of x2,
% and vice versa."

% We want to test the hypothesis that
% differences between real test and generated time series are larger than
% differences between pairs of real time series from the same subject.
% Then x1 is test vs sim, x2 is intra-subject.
% We want to test whether x1 tends to be larger than x2,
% so we test wither the CDF of x1 is smaller than that of x2.
[~, p_data_vs_generated_ps_gt_intra_ps] = kstest2( ps_rmse, same_subject_ps_rmse(:), 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|real-vs-generated and intra-subject PS RMSEs are from same distribution) = %g\n', ...
    p_data_vs_generated_ps_gt_intra_ps)
[~, p_data_vs_generated_fc_gt_intra_fc] = kstest2( fc_rmse(:), same_subject_fc_rmse(:), 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|real-vs-generated and intra-subject FC RMSEs are from same distribution) = %g\n', ...
    p_data_vs_generated_fc_gt_intra_fc)

% As a control, we also want to show that
% inter-subject differences tend to be greater than intra-subject ones.
% In this case x1 is inter-subject, and x2 is still intra-subject.
[~, p_inter_ps_gt_intra_ps] = kstest2( different_subject_ps_rmse(:), same_subject_ps_rmse(:), 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|inter-subject and intra-subject PS RMSEs are from same distribution) = %g\n', ...
    p_inter_ps_gt_intra_ps)
[~, p_inter_fc_gt_intra_fc] = kstest2( different_subject_fc_rmse(:), same_subject_fc_rmse(:), 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|inter-subject and intra-subject FC RMSEs are from same distribution) = %g\n', ...
    p_inter_fc_gt_intra_fc)

% Make boxplots to compare the RMSEs.

rmse_ts_fig = figure;
boxplot(ts_rmse)
ylabel('RMSE between real and simulated time series')
saveas(rmse_ts_fig, [figures_dir group '_single_ts_linear_model_time_series_rmse_box_plot_' settings_string '.fig'])

rmse_ps_fig = figure;
all_ps_rmses = [ 
    same_subject_ps_rmse(:);
    ps_rmse(:);
    different_subject_ps_rmse(:) 
    ];
all_ps_rmse_labels = [ 
    1*ones( numel(same_subject_ps_rmse), 1 ); 
    2*ones( numel(ps_rmse), 1 );
    3*ones( numel(different_subject_ps_rmse), 1 ) 
    ];
boxplot( log10(all_ps_rmses), all_ps_rmse_labels, 'Labels', {'intra-subject pairs', 'data vs generated', 'inter-subject pairs'} )
ylabel('log_{10}(RMSE between power spectra)')
saveas(rmse_ps_fig, [figures_dir group '_single_ts_linear_model_power_spectrum_rmse_box_plot_' settings_string '.fig'])

rmse_fc_fig = figure;
all_fc_rmses = [ 
    same_subject_fc_rmse(:);
    fc_rmse(:);
    different_subject_fc_rmse(:) 
    ];
all_fc_rmse_labels = [ 
    1*ones( numel(same_subject_fc_rmse), 1 ); 
    2*ones( numel(fc_rmse), 1 );
    3*ones( numel(different_subject_fc_rmse), 1 ) 
    ];
boxplot( all_fc_rmses, all_fc_rmse_labels, 'Labels', {'intra-subject pairs', 'data vs generated', 'inter-subject pairs'} )
ylabel('RMSE between functional connectivities above the diagonal')
saveas(rmse_ps_fig, [figures_dir group '_single_ts_linear_model_functional_connectivity_rmse_box_plot_' settings_string '.fig'])

% Plot the time series, spectra, and FCs of the last subject.

area_to_plot = 1;
area_str = sprintf('area_%u', area_to_plot);

% ts_ylim = [0.0 1.1];
% ts_ylim = [-200 250];
ts_ylim = [1.1*rand_min 1.1*rand_max];
ts_plot_fig = figure;
plot( time_series_times, data_ts(area_to_plot,:), '-r',...
    time_series_times, sim_ts(area_to_plot,:), '--g' )
legend({'data', 'sim with i.c. from data'})
ylim(ts_ylim)
xlabel('time (seconds)')
ylabel('BOLD signal')
saveas(ts_plot_fig, [figures_dir group '_example_generated_time_series_with_single_ts_linear_model_' settings_string '_' last_subject_id_str '_' area_str '.fig'])

ps_ylim = [-120 40];
% ps_ylim = [-50 65];
ps_plot_fig = figure;
plot( power_spectrum_frequencies,pow2db( get_power_spectrum(data_ts(area_to_plot,:), sampling_frequency)), '-r', ...
    power_spectrum_frequencies,pow2db( get_power_spectrum(sim_ts(area_to_plot,:), sampling_frequency)), '--g' )
legend({'data', 'sim with i.c. from data'})
ylim(ps_ylim)
xlabel('frequency (Hz)')
ylabel('power/frequency (decibels/Hz)')
saveas(ps_plot_fig, [figures_dir group '_example_generated_power_spectrum_with_single_ts_linear_model_' settings_string '_' last_subject_id_str '_' area_str '.fig'])

fc_plot_fig = figure;
subplot(1,2,1)
imshow( get_functional_connectivity(data_ts) )
title('data')
subplot(1,2,2)
imshow( get_functional_connectivity(sim_ts) )
title('sim with i.c. from data')
saveas(fc_plot_fig, [figures_dir group '_example_generated_functional_connectivity_single_ts_linear_model_' settings_string '_' last_subject_id_str '_' area_str '.fig'])
