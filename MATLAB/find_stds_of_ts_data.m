% MAKE_SC_WEIGHT_MASK
% Adam Craig, 2023-07-10
% Make a 360x360 weight matrix such that
% any place on the SC matrix that is always 0 in the training data is 0. 
% I am intentionally putting the variable num_sc_training_subjects
% in multiple places in this code,
% instead of just having num_subjects = num_sc_training_subjects,
% because we should not be using any other subject list to calculate this.
% Only use the training data.

hcp_data_header

all_ts = NaN(num_brain_areas, num_time_points, num_sc_training_subjects, time_series_per_subject);
for s = 1:num_sc_training_subjects
    subject_id = sc_training_subject_ids(s);
    for t = 1:time_series_per_subject
        ts_suffix = time_series_strings{t};
        ts_file = get_time_series_file(subject_id, ts_suffix);
        all_ts(:,:,s,t) = load_time_series(ts_file);
    end
end

ts_mean = mean(all_ts, 'all');
ts_std = std( all_ts, 0, 'all' );
fprintf('fMRI time series mean: %g, std dev: %g\n', ts_mean, ts_std)
