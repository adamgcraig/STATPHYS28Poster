% LEAVE_ONE_OUT_VALIDATE_LINEAR_MODEL_V2
% by Adam Craig, 2023-04-13.
% Demonstrate leave-one-out cross-validation using
% spectral power RMSE and functional connectivity RMSE
% with a simple linear model.
% In this version, we use min-max normalization on the data
% so that we can compare the results
% to a version that includes a nonlinear function with a range of [-1,+1].
% We also compute the correlations between
% weight matrices, FC matrices and SC matrices.

hcp_data_header

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

nonlinearity = @(v) v;
nonlinearity_inverse = @(v) v;
nonlinearity_name = 'detrend';
rescale_name = 'detrend';
rescale_fun = @(ts) detrend_time_series(ts);
rand_min = -10.0;
rand_max = 10.0;

% nonlinearity = @(v) v;
% nonlinearity_inverse = @(v) v;
% nonlinearity_name = 'std_mean_norm';
% rescale_name = 'std_mean_norm';
% rescale_fun = @(ts) std_mean_normalize_ts(ts);
% rand_min = -2.0;
% rand_max = 2.0;

% nonlinearity = @(v) v;
% nonlinearity_inverse = @(v) v;
% nonlinearity_name = 'min_max_norm';
% rescale_name = 'min_max_norm';
% rescale_fun = @(ts) rescale_ts( ts, -1.0, 1.0 );
% rand_min = -1.0;
% rand_max = 1.0;

% nonlinearity = @(v) v;
% nonlinearity_inverse = @(v) v;
% nonlinearity_name = 'min_max_norm_positive';
% rescale_name = 'min_max_norm_positive';
% rescale_fun = @(ts) rescale_ts( ts, 0.0, 1.0 );
% rand_min = 0.0;
% rand_max = 1.0;

use_sc = false;
use_sc_string = sprintf('use_sc_%u', use_sc);

print_status_update_if_time('starting code for leave-one-out cross-validation of linear model...')
% For this, always use the training subjects
% so that we still have the validation and testing subjects for later use.
time_series_set = cell(time_series_per_subject,1);
ps_set = cell(time_series_per_subject,1);
fc_set = cell(time_series_per_subject,1);
ts_nans = NaN(num_training_subjects, time_series_per_subject);
sc_vs_train_fc_corr = ts_nans;
W_vs_sc_corr = ts_nans;
W_vs_train_fc_corr = ts_nans;
train_vs_test_ps_rmse = ts_nans;
train_vs_test_fc_rmse = ts_nans;
train_vs_sim_train_ps_rmse = ts_nans;
train_vs_sim_train_fc_rmse = ts_nans;
test_vs_sim_test_ps_rmse = ts_nans;
test_vs_sim_test_fc_rmse = ts_nans;
test_vs_sim_rand_ps_rmse = ts_nans;
test_vs_sim_rand_fc_rmse = ts_nans;
has_sc = false(num_training_subjects,1);
for subject_index = 1:num_training_subjects
    subject_id = training_subject_ids(subject_index);
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
        time_series_file = get_time_series_file(subject_id, time_series_strings{time_series_index});
        data_ts = rescale_fun( load_time_series(time_series_file) );
        time_series_set{time_series_index} = data_ts;
        ps_set{time_series_index} = get_power_spectra_all_areas(data_ts, sampling_frequency);
        fc_set{time_series_index} = corr(data_ts');
        print_status_update_if_time( sprintf('loading subject %u of %u, time series %u of %u', ...
            subject_index, num_training_subjects, time_series_index, time_series_per_subject) )
    end
    for time_series_index = 1:time_series_per_subject
        % Train the model on 3 of the 4 time series.
        train_indices = (1:time_series_per_subject) ~= time_series_index;
        training_time_series_set = time_series_set(train_indices);
        test_ts = time_series_set{time_series_index};
        W = train_linear_model_on_time_series(training_time_series_set, use_y_intercept, sc_input, nonlinearity_inverse);
        abs_W = abs( W(:,1+use_y_intercept:end) );
        % Generate time series using the real time series initial state.
        sim_training_time_series_set = cell(time_series_per_subject-1,1);
        sim_training_ps_set = cell(time_series_per_subject-1,1);
        for training_ts_index = 1:time_series_per_subject-1
            training_ts = training_time_series_set{training_ts_index};
            sim_training_ts = generate_time_series_with_linear_model( W, training_ts(:,1), num_time_points, use_y_intercept, sc_input, nonlinearity );
            sim_training_time_series_set{training_ts_index} = sim_training_ts;
            sim_training_ps_set{training_ts_index} = get_power_spectra_all_areas(sim_training_ts, sampling_frequency);
        end
        sim_test_ts = generate_time_series_with_linear_model( W, test_ts(:,1), num_time_points, use_y_intercept, sc_input, nonlinearity );
        % Generate a time series with random initial conditions.
        rand_x0 = (rand_max-rand_min)*rand(num_brain_areas,1) + rand_min;
        sim_rand_ts = generate_time_series_with_linear_model(W, rand_x0, num_time_points, use_y_intercept, sc_input, nonlinearity);
        % Compute power spectra and functional connectivities.
        train_ps = mean( cat(3, ps_set{train_indices}), 3 );
        train_fc = get_functional_connectivity(  horzcat(training_time_series_set{:}) );
        test_ps = get_power_spectra_all_areas(test_ts, sampling_frequency);
        test_fc = get_functional_connectivity(test_ts);
        sim_train_ps = mean( cat(3, sim_training_ps_set{:}), 3 );
        sim_train_fc = corr( horzcat(sim_training_time_series_set{:})' );
        sim_test_ps = get_power_spectra_all_areas(sim_test_ts, sampling_frequency);
        sim_test_fc = get_functional_connectivity(sim_test_ts);
        sim_rand_ps = get_power_spectra_all_areas(sim_rand_ts, sampling_frequency);
        sim_rand_fc = get_functional_connectivity(sim_rand_ts);
        % Compute correlations.
        if current_subject_has_sc
            W_vs_sc_corr(subject_index, time_series_index) = get_fc_corr(abs_W, data_sc);
            sc_vs_train_fc_corr(subject_index, time_series_index) = get_fc_corr(data_sc, train_fc);
        end
        W_vs_train_fc_corr(subject_index, time_series_index) = get_fc_corr( abs_W, abs(train_fc) );
        % Compute RMSEs.
        % The power spectra matrices have no inherent symmetry or identity.
        % Compare the full spectra.
        train_vs_test_ps_rmse(subject_index, time_series_index) = get_rmse(train_ps, test_ps);
        train_vs_sim_train_ps_rmse(subject_index, time_series_index) = get_rmse(train_ps, sim_train_ps);
        test_vs_sim_test_ps_rmse(subject_index, time_series_index) = get_rmse(test_ps, sim_test_ps);
        test_vs_sim_rand_ps_rmse(subject_index, time_series_index) = get_rmse(test_ps, sim_rand_ps);
        % The FC matrix is symmetric with unit diagonal.
        % The RMSE between elements above the diagonal is more meaningful.
        train_vs_test_fc_rmse(subject_index, time_series_index) = get_upper_triangular_rmse(train_fc, test_fc);
        train_vs_sim_train_fc_rmse(subject_index, time_series_index) = get_upper_triangular_rmse(train_fc, sim_train_fc);
        test_vs_sim_test_fc_rmse(subject_index, time_series_index) = get_upper_triangular_rmse(test_fc, sim_test_fc);
        test_vs_sim_rand_fc_rmse(subject_index, time_series_index) = get_upper_triangular_rmse(test_fc, sim_rand_fc);
        % Print an update if enough time has passed.
        print_status_update_if_time( sprintf('testing against subject %u of %u, time series %u of %u', ...
            subject_index, num_training_subjects, time_series_index, time_series_per_subject) )
    end
end
last_subject_id_str = sprintf('subject_%u', subject_id);

settings_string = [ use_y_intercept_string '_' use_sc_string '_' nonlinearity_name ];

% Combine the stored results into a single table.
subject_id = reshape( repmat(training_subject_ids, 1, time_series_per_subject), [], 1 );
time_series = reshape( repmat(time_series_strings, num_training_subjects, 1), [], 1 );
has_sc = reshape( repmat(has_sc, 1, time_series_per_subject), [], 1 );
sc_vs_train_fc_corr = sc_vs_train_fc_corr(:);
W_vs_sc_corr = W_vs_sc_corr(:);
W_vs_train_fc_corr = W_vs_train_fc_corr(:);
train_vs_test_ps_rmse = train_vs_test_ps_rmse(:);
train_vs_test_fc_rmse = train_vs_test_fc_rmse(:);
train_vs_sim_train_ps_rmse = train_vs_sim_train_ps_rmse(:);
train_vs_sim_train_fc_rmse = train_vs_sim_train_fc_rmse(:);
test_vs_sim_test_ps_rmse = test_vs_sim_test_ps_rmse(:);
test_vs_sim_test_fc_rmse = test_vs_sim_test_fc_rmse(:);
test_vs_sim_rand_ps_rmse = test_vs_sim_rand_ps_rmse(:);
test_vs_sim_rand_fc_rmse = test_vs_sim_rand_fc_rmse(:);
results_table = table(subject_id, time_series, has_sc, ...
    sc_vs_train_fc_corr, W_vs_sc_corr, W_vs_train_fc_corr, ...
    train_vs_test_ps_rmse, train_vs_sim_train_ps_rmse, test_vs_sim_test_ps_rmse, test_vs_sim_rand_ps_rmse, ...
    train_vs_test_fc_rmse, train_vs_sim_train_fc_rmse, test_vs_sim_test_fc_rmse, test_vs_sim_rand_fc_rmse);
writetable(results_table, ['linear_model_training_leave_1_out_cross_validation_rmses_and_correlations_' settings_string '.csv'])

% Compare the FC and PS RMSE distributions to
% the distribution of RMSEs between FCs/PSs
% from different time series of the same subject.
% For comparison, compare the distribution of intra-subject RMSEs to
% the distribution of inter-subject RMSEs.

same_subject_ps_rmse = readmatrix([rmse_matrices_dir 'training_fmri_power_spectrum_rmses_same_subject_' rescale_name '.csv']);
different_subject_ps_rmse = readmatrix([rmse_matrices_dir 'training_fmri_power_spectrum_rmses_inter_subject_' rescale_name '.csv']);
same_subject_fc_rmse = readmatrix([rmse_matrices_dir 'training_fmri_functional_connectivity_rmses_same_subject_' rescale_name '.csv']);
different_subject_fc_rmse = readmatrix([rmse_matrices_dir 'training_fmri_functional_connectivity_rmses_inter_subject_' rescale_name '.csv']);

% Do statistical tests to
% see which RMSE distributions are different from each other.
% From https://ww2.mathworks.cn/help/stats/kstest2.html#namevaluepairarguments
% 'Tail', 'smaller'
% "Test the alternative hypothesis that the empirical cdf of x1 is smaller than the empirical cdf of x2."
% "If the data values in x1 tend to be larger than those in x2,
% the empirical distribution function of x1 tends to be smaller than that of x2,
% and vice versa."

% We want to test the hypothesis that
% differences between real test ts and generated time series
% are larger than
% differences between real training ts and generated time series
% Then x1 is test vs sim i.c., x2 is train vs sim
% We want to test whether x1 tends to be larger than x2,
% so we test wither the CDF of x1 is smaller than that of x2.
[~, p_test_vs_sim_ps_gt_train_vs_sim_fc] = kstest2( test_vs_sim_test_ps_rmse, train_vs_sim_train_ps_rmse, 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|test-vs-generated and train-vs-generated PS RMSEs are from same distribution) = %g\n', ...
    p_test_vs_sim_ps_gt_train_vs_sim_fc)
[~, p_test_vs_sim_fc_gt_train_vs_sim_fc] = kstest2( test_vs_sim_test_fc_rmse, train_vs_sim_train_fc_rmse, 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|test-vs-generated and train-vs-generated FC RMSEs are from same distribution) = %g\n', ...
    p_test_vs_sim_fc_gt_train_vs_sim_fc)

% We want to test the hypothesis that
% differences between real test ts and generated time series with random initial conditions
% are larger than
% differences between real test ts and generated time series with real initial conditions
% Then x1 is test vs sim with rand i.c., x2 is test vs sim with real i.c.
% We want to test whether x1 tends to be larger than x2,
% so we test wither the CDF of x1 is smaller than that of x2.
[~, p_test_vs_sim_rand_pc_gt_train_vs_sim_test_pc] = kstest2( test_vs_sim_rand_ps_rmse, test_vs_sim_test_ps_rmse, 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|test-vs-generated with rand i.c. and test-vs-generated with real i.c. PS RMSEs are from same distribution) = %g\n', ...
    p_test_vs_sim_rand_pc_gt_train_vs_sim_test_pc)
[~, p_test_vs_sim_rand_fc_gt_train_vs_sim_test_fc] = kstest2( test_vs_sim_rand_fc_rmse, test_vs_sim_test_fc_rmse, 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|test-vs-generated with rand i.c. and test-vs-generated with real i.c. FC RMSEs are from same distribution) = %g\n', ...
    p_test_vs_sim_rand_fc_gt_train_vs_sim_test_fc)

% We want to test the hypothesis that
% differences between real test and generated time series are larger than
% differences between pairs of real time series from the same subject.
% Then x1 is test vs sim, x2 is intra-subject.
% We want to test whether x1 tends to be larger than x2,
% so we test wither the CDF of x1 is smaller than that of x2.
[~, p_test_vs_generated_ps_gt_intra_ps] = kstest2( test_vs_sim_test_ps_rmse, same_subject_ps_rmse(:), 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|test-vs-generated and intra-subject PS RMSEs are from same distribution) = %g\n', ...
    p_test_vs_generated_ps_gt_intra_ps)
[~, p_test_vs_generated_fc_gt_intra_fc] = kstest2( test_vs_sim_test_fc_rmse(:), same_subject_fc_rmse(:), 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|test-vs-generated and intra-subject FC RMSEs are from same distribution) = %g\n', ...
    p_test_vs_generated_fc_gt_intra_fc)

% As a control, we also want to show that
% inter-subject differences tend to be greater than intra-subject ones.
% In this case x1 is inter-subject, and x2 is still intra-subject.
[~, p_inter_ps_gt_intra_ps] = kstest2( different_subject_ps_rmse(:), same_subject_ps_rmse(:), 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|inter-subject and intra-subject PS RMSEs are from same distribution) = %g\n', ...
    p_inter_ps_gt_intra_ps)
[~, p_inter_fc_gt_intra_fc] = kstest2( different_subject_fc_rmse(:), same_subject_fc_rmse(:), 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|inter-subject and intra-subject FC RMSEs are from same distribution) = %g\n', ...
    p_inter_fc_gt_intra_fc)

% Make a boxplot to compare the correlations.
corr_fig = figure;
boxplot( results_table{results_table.has_sc,{'sc_vs_train_fc_corr', 'W_vs_sc_corr', 'W_vs_train_fc_corr'}}, 'Labels', {'SC vs training data FC', 'W vs SC', 'W vs training data FC'} )
ylabel('Pearson correlation between absolute values of elements above diagonal')
saveas(corr_fig, [figures_dir 'linear_model_cross_validation_W_SC_FC_correlation_box_plot_' settings_string '.fig'])

% Make boxplots to compare the RMSEs.

rmse_ps_fig = figure;
all_ps_rmses = [ 
    same_subject_ps_rmse(:);
    train_vs_test_ps_rmse(:);
    train_vs_sim_train_ps_rmse(:);
    test_vs_sim_test_ps_rmse(:);
    test_vs_sim_rand_ps_rmse(:);
    different_subject_ps_rmse(:) 
    ];
all_ps_rmse_labels = [ 
    1*ones( numel(same_subject_ps_rmse), 1 ); 
    2*ones( numel(train_vs_test_ps_rmse), 1 );
    3*ones( numel(train_vs_sim_train_ps_rmse), 1 ); 
    4*ones( numel(test_vs_sim_test_ps_rmse), 1 ); 
    5*ones( numel(test_vs_sim_rand_ps_rmse), 1 ); 
    6*ones( numel(different_subject_ps_rmse), 1 ) 
    ];
boxplot( log10(all_ps_rmses), all_ps_rmse_labels, 'Labels', {'intra-subject pairs', 'mean train vs test', 'mean train vs mean sim train', 'test vs sim test', 'test vs sim rand i.c.', 'inter-subject pairs'} )
ylabel('log_{10}(RMSE between power spectra)')
saveas(rmse_ps_fig, [figures_dir 'linear_model_cross_validation_power_spectrum_rmse_box_plot_' settings_string '.fig'])

rmse_fc_fig = figure;
all_fc_rmses = [ 
    same_subject_fc_rmse(:);
    train_vs_test_fc_rmse(:);
    train_vs_sim_train_fc_rmse(:);
    test_vs_sim_test_fc_rmse(:);
    test_vs_sim_rand_fc_rmse(:);
    different_subject_fc_rmse(:) 
    ];
all_fc_rmse_labels = [ 
    1*ones( numel(same_subject_fc_rmse), 1 ); 
    2*ones( numel(train_vs_test_fc_rmse), 1 );
    3*ones( numel(train_vs_sim_train_fc_rmse), 1 ); 
    4*ones( numel(test_vs_sim_test_fc_rmse), 1 ); 
    5*ones( numel(test_vs_sim_rand_fc_rmse), 1 ); 
    6*ones( numel(different_subject_fc_rmse), 1 ) 
    ];
boxplot( all_fc_rmses, all_fc_rmse_labels, 'Labels', {'intra-subject pairs', 'mean train vs test', 'mean train vs mean sim train', 'test vs sim test', 'test vs sim rand i.c.', 'inter-subject pairs'} )
ylabel('RMSE between functional connectivities above the diagonal')
saveas(rmse_ps_fig, [figures_dir 'linear_model_cross_validation_functional_connectivity_rmse_box_plot_' settings_string '.fig'])

% Plot the time series, spectra, and FCs of the last subject.

area_to_plot = 1;
area_str = sprintf('area_%u', area_to_plot);

training_ts_1 = training_time_series_set{1};
training_ts_2 = training_time_series_set{2};
training_ts_3 = training_time_series_set{3};

sim_training_ts_1 = sim_training_time_series_set{1};
sim_training_ts_2 = sim_training_time_series_set{2};
sim_training_ts_3 = sim_training_time_series_set{3};

% ts_ylim = [0.0 1.1];
ts_ylim = [-200 250];
ts_plot_fig = figure;
subplot(2,2,1)
plot( time_series_times, training_ts_1(area_to_plot,:), '-r',...
    time_series_times, sim_training_ts_1(area_to_plot,:), '--g' )
legend({'data', 'sim with i.c. from data'})
ylim(ts_ylim)
xlabel('time (seconds)')
ylabel('BOLD signal')
title('training sequence 1')
subplot(2,2,2)
plot( time_series_times, training_ts_2(area_to_plot,:), '-r',...
    time_series_times, sim_training_ts_2(area_to_plot,:), '--g' )
legend({'data', 'sim with i.c. from data'})
ylim(ts_ylim)
xlabel('time (seconds)')
ylabel('BOLD signal')
title('training sequence 2')
subplot(2,2,3)
plot( time_series_times, training_ts_3(area_to_plot,:), '-r',...
    time_series_times, sim_training_ts_3(area_to_plot,:), '--g' )
legend({'data', 'sim with i.c. from data'})
ylim(ts_ylim)
xlabel('time (seconds)')
ylabel('BOLD signal')
title('training sequence 3')
subplot(2,2,4)
plot( time_series_times, test_ts(area_to_plot,:), '-r', ...
    time_series_times, sim_test_ts(area_to_plot,:), '--g', ...
    time_series_times, sim_rand_ts(area_to_plot,:), ':b' )
legend({'data', 'sim with i.c. from data', 'sim with random i.c.'})
ylim(ts_ylim)
xlabel('time (seconds)')
ylabel('BOLD signal')
title('testing sequence')
saveas(ts_plot_fig, [figures_dir 'example_3_vs_1_vs_generated_time_series_with_linear_model_' settings_string '_' last_subject_id_str '_' area_str '.fig'])

% ps_ylim = [-120 40];
ps_ylim = [-50 65];
ps_plot_fig = figure;
subplot(2,2,1)
plot( power_spectrum_frequencies,pow2db( get_power_spectrum(training_ts_1(area_to_plot,:), sampling_frequency)), '-r', ...
    power_spectrum_frequencies,pow2db( get_power_spectrum(sim_training_ts_1(area_to_plot,:), sampling_frequency)), '--g' )
legend({'data', 'sim with i.c. from data'})
ylim(ps_ylim)
xlabel('frequency (Hz)')
ylabel('power/frequency (decibels/Hz)')
title('training sequence 1')
subplot(2,2,2)
plot( power_spectrum_frequencies, pow2db(get_power_spectrum(training_ts_2(area_to_plot,:), sampling_frequency)), '-r', ...
    power_spectrum_frequencies,pow2db( get_power_spectrum(sim_training_ts_2(area_to_plot,:), sampling_frequency)), '--g' )
legend({'data', 'sim with i.c. from data'})
ylim(ps_ylim)
xlabel('frequency (Hz)')
ylabel('power/frequency (decibels/Hz)')
title('training sequence 2')
subplot(2,2,3)
plot( power_spectrum_frequencies, pow2db(get_power_spectrum(training_ts_3(area_to_plot,:), sampling_frequency)), '-r', ...
    power_spectrum_frequencies,pow2db( get_power_spectrum(sim_training_ts_3(area_to_plot,:), sampling_frequency)), '--g' )
legend({'data', 'sim with i.c. from data'})
ylim(ps_ylim)
xlabel('frequency (Hz)')
ylabel('power/frequency (decibels/Hz)')
title('training sequence 3')
subplot(2,2,4)
plot( power_spectrum_frequencies, pow2db(test_ps(area_to_plot,:)), '-r', ...
    power_spectrum_frequencies, pow2db(sim_test_ps(area_to_plot,:)), '--g', ...
    power_spectrum_frequencies, pow2db(sim_rand_ps(area_to_plot,:)), ':b' )
legend({'data', 'sim with i.c. from data', 'sim with random i.c.'})
ylim(ps_ylim)
xlabel('frequency (Hz)')
ylabel('power/frequency (decibels/Hz)')
title('testing sequence')
saveas(ps_plot_fig, [figures_dir 'example_3_vs_1_vs_generated_power_spectra_with_linear_model_' settings_string '_' last_subject_id_str '_' area_str '.fig'])

fc_plot_fig = figure;
subplot(2,5,1)
imshow( get_functional_connectivity(training_ts_1) )
title('training time series 1')
subplot(2,5,2)
imshow( get_functional_connectivity(training_ts_2) )
title('training time series 2')
subplot(2,5,3)
imshow( get_functional_connectivity(training_ts_3) )
title('training time series 3')
subplot(2,5,4)
imshow( test_fc )
title('testing time series')
subplot(2,5,6)
imshow( get_functional_connectivity(sim_training_ts_1) )
title('sim training time series 1')
subplot(2,5,7)
imshow( get_functional_connectivity(sim_training_ts_2) )
title('sim training time series 2')
subplot(2,5,8)
imshow( get_functional_connectivity(sim_training_ts_3) )
title('sim training time series 3')
subplot(2,5,9)
imshow( sim_test_fc )
title('sim testing time series')
subplot(2,5,10)
imshow( sim_rand_fc )
title('sim random i.c. time series')
saveas(fc_plot_fig, [figures_dir 'example_3_vs_1_vs_generated_functional_connectivity_with_linear_model_' settings_string '_' last_subject_id_str '_' area_str '.fig'])
