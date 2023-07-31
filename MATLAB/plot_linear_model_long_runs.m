% PLOT_LINEAR_MODEL_LONG_RUNS
% Adam Craig, 2023-06-08
% Check what happens when we run a linear model past 1200 time steps.

hcp_data_header

print_status_update_if_time('starting code for plotting long runs of linear models...')

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

use_sc = false;
use_sc_string = sprintf('use_sc_%u', use_sc);

extended_num_time_points = 2*num_time_points;
num_tp_string = sprintf('T_%u', extended_num_time_points);
settings_string = [ use_y_intercept_string '_' use_sc_string '_' nonlinearity_name '_' num_tp_string ];

group = 'training';
subject_ids = training_subject_ids;
% group = 'validation';
% subject_ids = validation_subject_ids;
% group = 'testing';
% subject_ids = testing_subject_ids;

num_subjects = numel(subject_ids);

subject_index = 50;% num_subjects;

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
training_time_series = cell(time_series_per_subject,1);
training_ps = cell(time_series_per_subject,1);
training_fc = cell(time_series_per_subject,1);
for time_series_index = 1:time_series_per_subject
    ts_suffix = time_series_strings{time_series_index};
    time_series_file = get_time_series_file(subject_id, ts_suffix);
    training_time_series{time_series_index} = rescale_fun( load_time_series(time_series_file) );
    training_ps{time_series_index} = get_power_spectra_all_areas(data_ts, sampling_frequency);
    training_fc{time_series_index} = get_functional_connectivity(data_ts);
end
extended_time_series_times = 0:time_series_resolution:time_series_resolution*(extended_num_time_points-1);
area_to_plot = 1;
area_str = sprintf('area_%u', area_to_plot);
for num_training_ts = 1:time_series_per_subject
    training_ts_string_plot = strjoin( time_series_plot_text_strings(1:num_training_ts), ', ' );
    training_ts_string_file_name = strjoin( time_series_strings(1:num_training_ts), '_' );
    W = train_linear_model_on_time_series( training_time_series(1:num_training_ts), use_y_intercept, sc_input, nonlinearity_inverse );
    for time_series_index = 1:time_series_per_subject
        ts_suffix = time_series_strings{time_series_index};
        ts_plot_string = time_series_plot_text_strings{time_series_index};
        ts_string = sprintf('%u_%s', subject_id, ts_suffix);
        data_ts = training_time_series{time_series_index};
        data_ps = training_ps{time_series_index};
        data_fc = training_fc{time_series_index};
        sim_ts = generate_time_series_with_linear_model( W, data_ts(:,1), extended_num_time_points, use_y_intercept, sc_input, nonlinearity );
        % sim_ps = get_power_spectra_all_areas(sim_ts, sampling_frequency);
        % sim_fc = get_functional_connectivity(sim_ts);
        % ts_rmse = get_rmse( sim_ts(:,2:end), data_ts(:,2:end) );
        % ps_rmse = get_rmse(sim_ps, data_ps);
        % fc_rmse = get_upper_triangular_rmse(sim_fc, data_fc);

        % ts_ylim = [0.0 1.1];
        % ts_ylim = [-200 250];
        ts_ylim = [1.1*rand_min 1.1*rand_max];
        ts_plot_fig = figure;
        plot( time_series_times, data_ts(area_to_plot,:), '-r',...
            extended_time_series_times, sim_ts(area_to_plot,:), '--g' )
        legend({'data', 'sim with i.c. from data'})
        ylim(ts_ylim)
        xlabel('time (seconds)')
        ylabel('BOLD signal')
        title( sprintf('model of subject %u, trained on %s tested on %s', subject_id, training_ts_string_plot, ts_plot_string) )
        fig_name = [figures_dir group '_example_extended_generated_time_series_' settings_string '_' ts_string '_trained_on_' training_ts_string_file_name '_' area_str];
        saveas(ts_plot_fig, [fig_name '.fig'])
        saveas(ts_plot_fig, [fig_name '.png'])

        % ps_ylim = [-120 40];
        % ps_ylim = [-50 65];
        % ps_plot_fig = figure;
        % plot( power_spectrum_frequencies,pow2db( get_power_spectrum(data_ts(area_to_plot,:), sampling_frequency)), '-r', ...
        %     power_spectrum_frequencies,pow2db( get_power_spectrum(sim_ts(area_to_plot,:), sampling_frequency)), '--g' )
        % legend({'data', 'sim with i.c. from data'})
        % ylim(ps_ylim)
        % xlabel('frequency (Hz)')
        % ylabel('power/frequency (decibels/Hz)')
        % saveas(ps_plot_fig, [figures_dir group '_example_generated_power_spectrum_with_single_ts_linear_model_' settings_string '_' last_subject_id_str '_' area_str '.fig'])

        % fc_plot_fig = figure;
        % subplot(1,2,1)
        % imshow( get_functional_connectivity(data_ts) )
        % title('data')
        % subplot(1,2,2)
        % imshow( get_functional_connectivity(sim_ts) )
        % title('sim with i.c. from data')
        % saveas(ts_plot_fig, [figures_dir group '_example_extended_generated_fc_' settings_string '_' ts_string '_trainded_on_' sprintf('%u', num_training_subjects) '.fig'])

    end
end
