% PLOT_FFT_COMPONENTS_VS_RMSE.M
% by Adam Craig, 2023-07-26
% Plot the relationship between
% number of Fourier coefficients used to reconstruct a HCP fMRI time series
% and
% root mean squared error (RMSE) between the original and reconstruction.
% Use the training data set.

% hcp_data_header

% area_to_plot = 1;
% subject_id = training_subject_ids(1);
% ts_string = time_series_strings{1};
% ts_file = get_time_series_file(subject_id, ts_string);
% ts = load_time_series(ts_file);

% ts_fft = fft(ts, num_time_points, 2);
% midpoint = num_time_points/2 + 1;
% ts_fft_zeroed = ts_fft;
% freqs_to_erase = 300;
% ts_fft_zeroed(:,midpoint-freqs_to_erase+1:midpoint+freqs_to_erase-1) = 0;
% freq_component = 1:num_time_points;
% figure
% plot(  freq_component, abs( ts_fft(area_to_plot,:) ), '-r', freq_component, abs( ts_fft_zeroed(area_to_plot,:) ), '--g'  )
% legend({'original', 'partially zeroed'})
% freq_component = 1:num_time_points;
% figure
% plot(  freq_component, abs( ts_fft(area_to_plot,:) ), '-r'  )
% xlabel('magnitude')
% ylabel('frequency component index')

% ts_recon = ifft(ts_fft, num_time_points, 2, 'symmetric');
% ts_rmse = get_rmse(ts, ts_recon);
% time_point = 1:num_time_points;
% figure
% plot( time_point, ts(area_to_plot,:), '-r', time_point, ts_recon(area_to_plot,:), '--g' )
% legend({'original', 'reconstruction'})
% fprintf( 'original RMSE: %g\n', ts_rmse )

subjects_ids = training_subject_ids;
num_subjects = numel(subjects_ids);
ts_per_subject = numel(time_series_strings);
num_freqs = floor(num_time_points/2);
freqs_to_erase_choices = 0:num_freqs;
num_choices = numel(freqs_to_erase_choices);
ts_rmse = zeros(num_choices, num_subjects, ts_per_subject);
midpoint = num_time_points/2 + 1;
for subject_index = 1:num_subjects
    subject_id = subjects_ids(subject_index);
    for ts_index = 1:ts_per_subject
        ts_string = time_series_strings{ts_index};
        ts_file = get_time_series_file(subject_id, ts_string);
        ts = min_max_normalize_ts( load_time_series(ts_file) );
        ts_fft = fft(ts, num_time_points, 2);
        for choice_index = 1:num_choices
            freqs_to_erase = freqs_to_erase_choices(choice_index);
            ts_fft(:,midpoint-freqs_to_erase+1:midpoint+freqs_to_erase-1) = 0;
            ts_recon = ifft(ts_fft, num_time_points, 2, 'symmetric');
            ts_rmse(choice_index,subject_index,ts_index) = get_rmse(ts, ts_recon);
        end
    end
end
% figure
% plot( time_point, ts(area_to_plot,:), '-r', time_point, ts_recon(area_to_plot,:), '--g' )
% legend({'original', 'reconstruction'})
% fprintf( 'finale RMSE: %g\n', ts_rmse(end) )

min_ts_rmse = min( min(ts_rmse, [], 3), [], 2 );
mean_ts_rmse = mean( mean(ts_rmse, 3), 2 );
max_ts_rmse = max( max(ts_rmse, [], 3), [], 2 );
ts_rmse_neg_diff = mean_ts_rmse-min_ts_rmse;
ts_rmse_pos_diff = max_ts_rmse-mean_ts_rmse;

% figure
% errorbar(freqs_to_erase_choices, mean_ts_rmse, ts_rmse_neg_diff, ts_rmse_pos_diff)
% xlabel('frequencies erased')
% ylabel('RMSE (range around mean)')

sparse_indices = 1:10:num_freqs-1;
figure
errorbar( freqs_to_erase_choices(sparse_indices), mean_ts_rmse(sparse_indices), ts_rmse_neg_diff(sparse_indices), ts_rmse_pos_diff(sparse_indices) )
xlabel('frequencies erased')
ylabel('RMSE (range around mean)')

sparse_indices = 1:10:500;% num_freqs-1;
figure
errorbar( freqs_to_erase_choices(sparse_indices), mean_ts_rmse(sparse_indices), ts_rmse_neg_diff(sparse_indices), ts_rmse_pos_diff(sparse_indices) )
xlabel('frequencies erased')
ylabel('RMSE (range around mean)')
