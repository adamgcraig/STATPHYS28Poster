% COMPUTE_ALL_DATA_POWER_SPECTRA_AND_PAIRWISE_RMSES
% by Adam Craig, 2023-04-13.
% Compute the power spectra of the HCP fMRI time series,
% and calculate the RMSE between each pair of time series.

print_status_update_if_time('starting code for compute_all_data_power_spectra_and_pairwise_rmses')

hcp_data_header

% rescale_name = 'raw';
% rescale_fun = @(ts) ts;

% rescale_name = 'zero_mean';
% rescale_fun = @(ts) ts - mean(ts,'all');

% rescale_name = 'detrend';
% rescale_fun = @(ts) detrend_time_series(ts);

% rescale_name = 'std_mean_norm';
% rescale_fun = @(ts) std_mean_normalize_ts(ts);

rescale_name = 'min_max_norm';
rescale_fun = @(ts) rescale_ts( ts, -1.0, 1.0 );

% rescale_name = 'min_max_norm_positive';
% rescale_fun = @(ts) rescale_ts( ts, 0.0, 1.0 );

power_spectra = NaN(num_brain_areas, num_power_spectrum_frequencies);
% Here we can swap in
% training_subject_ids, validation_subject_ids or testing_subject_ids.
% output_file_prefix = 'training';
% subject_ids = training_subject_ids;
% output_file_prefix = 'validation';
% subject_ids = validation_subject_ids;
output_file_prefix = 'testing';
subject_ids = testing_subject_ids;
num_subjects = numel(subject_ids);
power_spectrum_files = cell(num_subjects, time_series_per_subject);
disp('computing power spectra...')
for subject_index = 1:num_subjects
    subject_id = subject_ids(subject_index);
    for time_series_index = 1:time_series_per_subject
        time_series_string = time_series_strings{time_series_index};
        time_series_file = get_time_series_file(subject_id, time_series_string);
        time_series = rescale_fun( load_time_series(time_series_file) );
        for brain_area_index = 1:num_brain_areas
            power_spectra(brain_area_index,:) = get_power_spectrum( time_series(brain_area_index,:), sampling_frequency );
            print_status_update_if_time( sprintf('subject %u of %u, time series %u of %u, area %u of %u', ...
                subject_index, num_subjects, time_series_index, time_series_per_subject, brain_area_index, num_brain_areas) )
        end
        power_spectrum_file = get_power_spectrum_file(subject_id, time_series_string);
        write_binary_data(power_spectra, power_spectrum_file)
        power_spectrum_files{subject_index, time_series_index} = power_spectrum_file;
    end
end

disp('computing pairwise RMSEs between power spectra...')
num_time_series = numel(power_spectrum_files);
power_spectrum_rmses = zeros(num_time_series, num_time_series);
for time_series_1_index = 1:num_time_series-1
    power_spectrum_1 = load_data_from_binary(power_spectrum_files{time_series_1_index}, num_brain_areas, num_power_spectrum_frequencies);
    for time_series_2_index = time_series_1_index+1:num_time_series
        power_spectrum_2 = load_data_from_binary(power_spectrum_files{time_series_2_index}, num_brain_areas, num_power_spectrum_frequencies);
        rmse = get_rmse(power_spectrum_1, power_spectrum_2);
        power_spectrum_rmses(time_series_1_index, time_series_2_index) = rmse;
        power_spectrum_rmses(time_series_2_index, time_series_1_index) = rmse;
        print_status_update_if_time( sprintf('time series 1 %u, time series 2 %u', time_series_1_index, time_series_2_index) )
    end
end
writematrix(power_spectrum_rmses, [rmse_matrices_dir output_file_prefix '_fmri_power_spectrum_rmses_' rescale_name '.csv'])

disp('comparing same-subject RMSEs to inter-subject RMSEs...')
subject_id_for_time_series = repmat(subject_ids, 1, time_series_per_subject);
is_same_subject = reshape(subject_id_for_time_series, [], 1) == reshape(subject_id_for_time_series, 1, []);
% Select only values in the upper triangle above the diagonal.
% That way, we only include the RMSE between each pair of different spectra
% once.
ut_power_spectrum_rmses = get_upper_triangular_elements(power_spectrum_rmses);
ut_is_same_subject = get_upper_triangular_elements(is_same_subject);
same_subject_power_spectrum_rmses = ut_power_spectrum_rmses(ut_is_same_subject);
different_subject_power_spectrum_rmses = ut_power_spectrum_rmses(~ut_is_same_subject);
writematrix(same_subject_power_spectrum_rmses, [rmse_matrices_dir output_file_prefix '_fmri_power_spectrum_rmses_same_subject_' rescale_name '.csv'])
writematrix(different_subject_power_spectrum_rmses, [rmse_matrices_dir output_file_prefix '_fmri_power_spectrum_rmses_inter_subject_' rescale_name '.csv'])
[~, ttest_p_value] = ttest2(same_subject_power_spectrum_rmses, different_subject_power_spectrum_rmses);
fprintf('p-value of 2-sample t-test between same-subject and inter-subject power spectrum RMSEs: %g\n', ttest_p_value)

box_fig = figure;
boxplot(ut_power_spectrum_rmses, ut_is_same_subject, 'Labels', {'inter-subject', 'same-subject'})
ylabel('RMSE between fMRI power spectra of all brain areas')
saveas(box_fig, [figures_dir 'inter_vs_same_subject_power_spectrum_pairwise_rmse_boxplots_' rescale_name '.fig'])

hist_fig = figure;
hold on
histogram(different_subject_power_spectrum_rmses)
histogram(same_subject_power_spectrum_rmses)
legend({'inter-subject', 'same-subject'})
xlabel('RMSE between fMRI power spectra of all brain areas')
ylabel('number of pairs')
saveas(hist_fig, [figures_dir 'inter_vs_same_subject_power_spectrum_pairwise_rmse_histograms_' rescale_name '.fig'])

print_status_update_if_time( 'done' )
