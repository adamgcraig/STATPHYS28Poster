% LEAVE_ONE_OUT_VALIDATE_LINEAR_MODEL
% by Adam Craig, 2023-04-13.
% Demonstrate leave-one-out cross-validation using
% spectral power RMSE and functional connectivity RMSE
% with a simple linear model.

hcp_data_header

print_status_update_if_time('starting code for leave-one-out cross-validation of linear model...')
% For this, always use the training subjects
% so that we still have the validation and testing subjects for later use.
test_vs_sim_ps_rmses = NaN(num_training_subjects, time_series_per_subject);
test_vs_sim_fc_rmses = NaN(num_training_subjects, time_series_per_subject);
time_series_set = cell(time_series_per_subject,1);
for subject_index = 1:num_training_subjects
    subject_id = training_subject_ids(subject_index);
    for time_series_index = 1:time_series_per_subject
        time_series_file = get_time_series_file(subject_id, time_series_strings{time_series_index});
        time_series_set{time_series_index} = load_time_series(time_series_file);
        print_status_update_if_time( sprintf('loading subject %u of %u, time series %u of %u', ...
            subject_index, num_training_subjects, time_series_index, time_series_per_subject) )
    end
    for time_series_index = 1:time_series_per_subject
        W = train_linear_model_on_time_series(  time_series_set( (1:time_series_per_subject) ~= time_series_index )  );
        test_time_series = time_series_set{time_series_index};
        test_power_spectra = get_power_spectra_all_areas(test_time_series, sampling_frequency);
        test_fc = get_functional_connectivity(test_time_series);
        sim_time_series = generate_time_series_with_linear_model( W, test_time_series(:,1), num_time_points );
        sim_fc = get_functional_connectivity(sim_time_series);
        sim_power_spectra = get_power_spectra_all_areas(sim_time_series, sampling_frequency);
        % The power spectra matrices have no inherent symmetry or identity.
        % Compare the full spectra.
        test_vs_sim_ps_rmses(subject_index, time_series_index) = get_rmse(test_power_spectra, sim_power_spectra);
        % The FC matrix is symmetric with unit diagonal.
        % The RMSE between elements above the diagonal is more meaningful.
        test_vs_sim_fc_rmses(subject_index, time_series_index) = get_upper_triangular_rmse(test_fc, sim_fc);
        print_status_update_if_time( sprintf('testing against subject %u of %u, time series %u of %u', ...
            subject_index, num_training_subjects, time_series_index, time_series_per_subject) )
    end  
end
writematrix(test_vs_sim_ps_rmses, 'linear_model_training_leave_1_out_cross_validation_power_spectrum_rmses.csv')
writematrix(test_vs_sim_fc_rmses, 'linear_model_training_leave_1_out_cross_validation_functional_connectivity_rmses.csv')

same_subject_power_spectrum_rmses = readmatrix([rmse_matrices_dir 'training_fmri_power_spectrum_rmses_same_subject.csv']);
different_subject_power_spectrum_rmses = readmatrix([rmse_matrices_dir 'training_fmri_power_spectrum_rmses_inter_subject.csv']);
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
[~, p_test_vs_generated_ps_gt_intra_ps] = kstest2( test_vs_sim_ps_rmses(:), same_subject_power_spectrum_rmses(:), 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|test-vs-generated and intra-subject power spectrum RMSEs are from same distribution) = %g\n', ...
    p_test_vs_generated_ps_gt_intra_ps)
% As a control, we also want to show that
% inter-subject differences tend to be greater than intra-subject ones.
% In this case x1 is inter-subject, and x2 is still intra-subject.
[~, p_inter_ps_gt_intra_ps] = kstest2( different_subject_power_spectrum_rmses(:), same_subject_power_spectrum_rmses(:), 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|inter-subject and intra-subject power spectrum RMSEs are from same distribution) = %g\n', ...
    p_inter_ps_gt_intra_ps)

rmse_ps_fig = figure;
all_ps_rmses = [ 
    same_subject_power_spectrum_rmses(:); 
    test_vs_sim_ps_rmses(:); 
    different_subject_power_spectrum_rmses(:) 
    ];
all_ps_rmse_labels = [ 
    1*ones( numel(same_subject_power_spectrum_rmses), 1 ); 
    2*ones( numel(test_vs_sim_ps_rmses), 1 ); 
    3*ones( numel(different_subject_power_spectrum_rmses), 1 ) 
    ];
boxplot( log10(all_ps_rmses), all_ps_rmse_labels, 'Labels', {'intra-subject pairs', 'generated vs left-out', 'inter-subject pairs'} )
ylabel('log_{10}(RMSE between power spectra)')
saveas(rmse_ps_fig, [figures_dir 'linear_model_cross_validation_power_spectrum_rmse_box_plot.fig'])

% Same explanation above applies for testing functional connectivity RMSEs.
same_subject_functional_connectivity_rmses = readmatrix([rmse_matrices_dir 'training_fmri_functional_connectivity_rmses_same_subject.csv']);
different_subject_functional_connectivity_rmses = readmatrix([rmse_matrices_dir 'training_fmri_functional_connectivity_rmses_inter_subject.csv']);
[~, p_test_vs_generated_fc_gt_intra_fc] = kstest2( test_vs_sim_fc_rmses(:), same_subject_functional_connectivity_rmses(:), 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|test-vs-generated and intra-subject functional connectivity RMSEs are from same distribution) = %g\n', ...
    p_test_vs_generated_fc_gt_intra_fc)
[~, p_inter_fc_gt_intra_fc] = kstest2( different_subject_functional_connectivity_rmses(:), same_subject_functional_connectivity_rmses(:), 'Tail', 'smaller' );
fprintf('P(this or more extreme difference|inter-subject and intra-subject functional connectivity RMSEs are from same distribution) = %g\n', ...
    p_inter_fc_gt_intra_fc)

rmse_fc_fig = figure;
all_fc_rmses = [ 
    same_subject_functional_connectivity_rmses(:); 
    test_vs_sim_fc_rmses(:); 
    different_subject_functional_connectivity_rmses(:) 
    ];
all_fc_rmse_labels = [ 
    1*ones( numel(same_subject_functional_connectivity_rmses), 1 ); 
    2*ones( numel(test_vs_sim_fc_rmses), 1 ); 
    3*ones( numel(different_subject_functional_connectivity_rmses), 1 ) 
    ];
boxplot( all_fc_rmses, all_fc_rmse_labels, 'Labels', {'intra-subject pairs', 'generated vs left-out', 'inter-subject pairs'} )
ylabel('RMSE between functional connectivites')
saveas(rmse_fc_fig, [figures_dir 'linear_model_cross_validation_functional_connectivity_rmse_box_plot.fig'])

area_to_plot = 1;
training_ts_1 = time_series_set{1};
training_ts_2 = time_series_set{2};
training_ts_3 = time_series_set{3};

% Putting them all on the same axes gets too messy.
% ts_plot_fig = figure;
% plot( time_series_times, training_ts_1(area_to_plot,:), '-c', ...
%     time_series_times, training_ts_2(area_to_plot,:), '-b', ...
%     time_series_times, training_ts_3(area_to_plot,:), '-m', ...
%     time_series_times, test_time_series(area_to_plot,:), '-g', ...
%     time_series_times, sim_time_series(area_to_plot,:), '-r' )
% legend({'training 1', 'training 2', 'training 3', 'testing', 'generated'})
% xlabel('time (seconds)')
% ylabel('BOLD signal')
% saveas(ts_plot_fig, [figures_dir 'example_3_vs_1_vs_generated_time_series_with_linear_model.fig'])

ts_ylim = [-200 300];
ts_plot_fig = figure;
subplot(2,2,1)
plot( time_series_times, training_ts_1(area_to_plot,:), '-c' )
ylim(ts_ylim)
xlabel('time (seconds)')
ylabel('BOLD signal')
title('training sequence 1')
subplot(2,2,2)
plot( time_series_times, training_ts_2(area_to_plot,:), '-b' )
ylim(ts_ylim)
xlabel('time (seconds)')
ylabel('BOLD signal')
title('training sequence 2')
subplot(2,2,3)
plot( time_series_times, training_ts_3(area_to_plot,:), '-m' )
ylim(ts_ylim)
xlabel('time (seconds)')
ylabel('BOLD signal')
title('training sequence 3')
subplot(2,2,4)
plot( time_series_times, test_time_series(area_to_plot,:), '-g', ...
    time_series_times, sim_time_series(area_to_plot,:), '-r' )
legend({'testing', 'generated'})
ylim(ts_ylim)
xlabel('time (seconds)')
ylabel('BOLD signal')
title('generated sequence and testing sequence')
saveas(ts_plot_fig, [figures_dir 'example_3_vs_1_vs_generated_time_series_with_linear_model_separate.fig'])

% ps_plot_fig = figure;
% plot( power_spectrum_frequencies,pow2db( get_power_spectrum(training_ts_1(area_to_plot,:), sampling_frequency)), '-c', ...
%     power_spectrum_frequencies, pow2db(get_power_spectrum(training_ts_2(area_to_plot,:), sampling_frequency)), '-b', ...
%     power_spectrum_frequencies, pow2db(get_power_spectrum(training_ts_3(area_to_plot,:), sampling_frequency)), '-m', ...
%     power_spectrum_frequencies, pow2db(test_power_spectra(area_to_plot,:)), '-g', ...
%     power_spectrum_frequencies, pow2db(sim_power_spectra(area_to_plot,:)), '-r' )
% legend({'training 1', 'training 2', 'training 3', 'testing', 'generated'})
% xlabel('frequency (Hz)')
% ylabel('power/frequency (decibels/Hz)')
% saveas(ps_plot_fig, [figures_dir 'example_3_vs_1_vs_generated_power_spectra_with_linear_model.fig'])

ps_ylim = [-55 65];
ps_plot_fig = figure;
subplot(2,2,1)
plot( power_spectrum_frequencies,pow2db( get_power_spectrum(training_ts_1(area_to_plot,:), sampling_frequency)), '-c' )
ylim(ps_ylim)
xlabel('frequency (Hz)')
ylabel('power/frequency (decibels/Hz)')
title('training sequence 1')
subplot(2,2,2)
plot( power_spectrum_frequencies, pow2db(get_power_spectrum(training_ts_2(area_to_plot,:), sampling_frequency)), '-b' )
ylim(ps_ylim)
xlabel('frequency (Hz)')
ylabel('power/frequency (decibels/Hz)')
title('training sequence 2')
subplot(2,2,3)
plot( power_spectrum_frequencies, pow2db(get_power_spectrum(training_ts_3(area_to_plot,:), sampling_frequency)), '-m' )
ylim(ps_ylim)
xlabel('frequency (Hz)')
ylabel('power/frequency (decibels/Hz)')
title('training sequence 3')
subplot(2,2,4)
plot( power_spectrum_frequencies, pow2db(test_power_spectra(area_to_plot,:)), '-g', ...
    power_spectrum_frequencies, pow2db(sim_power_spectra(area_to_plot,:)), '-r' )
legend({'testing', 'generated'})
ylim(ps_ylim)
xlabel('frequency (Hz)')
ylabel('power/frequency (decibels/Hz)')
title('generated sequence vs testing sequence')
saveas(ps_plot_fig, [figures_dir 'example_3_vs_1_vs_generated_power_spectra_with_linear_model_separate.fig'])

fc_plot_fig = figure;
subplot(2,3,1)
imshow( get_functional_connectivity(training_ts_1) )
title('training time series 1')
subplot(2,3,2)
imshow( get_functional_connectivity(training_ts_2) )
title('training time series 2')
subplot(2,3,3)
imshow( get_functional_connectivity(training_ts_3) )
title('training time series 3')
subplot(2,3,4)
imshow( test_fc )
title('testing time series')
subplot(2,3,5)
imshow( sim_fc )
title('generated time series')
saveas(fc_plot_fig, [figures_dir 'example_3_vs_1_vs_generated_functional_connectivity_with_linear_model.fig'])

