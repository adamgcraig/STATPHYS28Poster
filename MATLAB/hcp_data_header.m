% HCP_DATA_HEADER
% by Adam Craig, 2023-04-13.
% Call this within other scripts to set environment variables we need
% in order to use the Human Connectome Project (HCP) data.

% Data directory. Update this to reflect your local file system.
data_dir = 'E:\HCP_data\';
% Optionally, we can send output to a separate directory.
output_dir = data_dir;

% Information about the brain regions.
% These come from the parcellation described in
% Glasser, M. F., Coalson, T. S., Robinson, E. C., Hacker, C. D., Harwell,
% J., Yacoub, E., ... & Van Essen, D. C. (2016).
% A multi-modal parcellation of human cerebral cortex.
% Nature, 536(7615), 171-178.

% MRI data are mapped onto a standard Atlas where the areas have fixed
% (x,y,z) coordinates and abreviations.
brain_area_table_file = [data_dir 'MMP360coordinator.xlsx'];
brain_area_table = get_brain_area_table(brain_area_table_file);
num_brain_areas = size(brain_area_table,1);

% Each subject has 1 area feature file.
% It is a 360 area by 4 feature matrix.
% We store it in a binary file of 64-bit floating point numbers.
brain_area_feature_names = {'thickness', 'myelination', 'curvature', 'sulcus depth'};
% The 4 features are
features_per_brain_area = numel(brain_area_feature_names);
brain_area_feature_dir =[data_dir 'anatomy_binaries\'];
get_brain_area_feature_file = @(subject_id) [brain_area_feature_dir sprintf('anatomy_%u.bin', subject_id)];

% Each subject has 1 structural connectivity file.
% It is a 360 area by 360 area symmetric, nonnegative matrix.
% We store it in a binary file of 64-bit floating point numbers.
structural_connectivity_dir = [data_dir 'dtMRI_binaries\'];
get_structural_connectivity_file = @(subject_id) [structural_connectivity_dir sprintf('sc_%u.bin', subject_id)];

% Each preprocessed fMRI time series is 360 area by 1200 time point matrix.
% We store it in a binary file of 64-bit floating point numbers.
% Each subject has 4 time series named with the suffixes below.
% The number indicates
% whether the scan was taken on the 1st or 2nd session.
% The letter pair indicates
% whether the scan was taken from left to right or from right to left.
time_series_strings = {'1_LR' '1_RL' '2_LR' '2_RL'};
% Use these versions in plot titles, legends, etc.
time_series_plot_text_strings = {'1LR' '1RL' '2LR' '2RL'};
time_series_per_subject = numel(time_series_strings);
time_series_dir = [data_dir 'fMRI_ts_binaries\'];
get_time_series_file = @(subject_id, time_series_string) [time_series_dir sprintf('ts_%u_%s.bin', subject_id, time_series_string)];
num_time_points = 1200;
time_series_resolution = 0.72;% seconds between samples
time_series_times = 0:time_series_resolution:time_series_resolution*(num_time_points-1);

% IDs of subjects for whom we have full time series
% and area features,
% shuffled and split into training, validation, and testing sets.
% There should be 838 subjects total.
training_subject_file = [data_dir 'training_subject_ids.txt'];
training_subject_ids = readmatrix(training_subject_file);
num_training_subjects = numel(training_subject_ids);
validation_subject_file = [data_dir 'validation_subject_ids.txt'];
validation_subject_ids = readmatrix(validation_subject_file);
num_validation_subjects = numel(validation_subject_ids);
testing_subject_file = [data_dir 'testing_subject_ids.txt'];
testing_subject_ids = readmatrix(testing_subject_file);
num_testing_subjects = numel(testing_subject_ids);

% IDs of subjects for whom we have full time series,
% area features, and structural connectivity,
% shuffled and split into training, validation, and testing sets.
% There should be 838 subjects total.
sc_training_subject_file = [data_dir 'sc_training_subject_ids.txt'];
sc_training_subject_ids = readmatrix(sc_training_subject_file);
num_sc_training_subjects = numel(sc_training_subject_ids);
sc_validation_subject_file = [data_dir 'sc_validation_subject_ids.txt'];
sc_validation_subject_ids = readmatrix(sc_validation_subject_file);
num_sc_validation_subjects = numel(sc_validation_subject_ids);
sc_testing_subject_file = [data_dir 'sc_testing_subject_ids.txt'];
sc_testing_subject_ids = readmatrix(sc_testing_subject_file);
num_sc_testing_subjects = numel(sc_testing_subject_ids);

% Set up some directories for different kinds of output.

% for power spectra
power_spectrum_dir = [output_dir 'power_spectra\'];
get_power_spectrum_file = @(subject_id, time_series_string) [power_spectrum_dir sprintf('spectrum_%u_%s.bin', subject_id, time_series_string)];
sampling_frequency = 1/time_series_resolution;
power_spectrum_frequencies = 0:sampling_frequency/num_time_points:sampling_frequency/2;
num_power_spectrum_frequencies = numel(power_spectrum_frequencies);

functional_connectivity_dir = [output_dir 'functional_connectivity\'];
get_functional_connectivity_file = @(subject_id, time_series_string) [functional_connectivity_dir sprintf('fc_%u_%s.bin', subject_id, time_series_string)];

model_dir = [output_dir 'models\'];
single_ts_linear_model_dir = [model_dir 'single_ts_linear\'];
generated_time_series_dir = [output_dir 'generated_time_series\'];
rmse_matrices_dir = [output_dir 'rmse_matrices\'];
figures_dir = [output_dir 'figures\'];

interpolate_ts = @(ts, interp_factor) spline( 1:size(ts,2), ts, (1/interp_factor):(1/interp_factor):size(ts,2) );
std_mean_normalize_ts = @(ts) ( ts - mean(ts,'all') )./std(ts,0,'all');
min_max_normalize_ts = @(ts) ( ts - min(ts,[],'all') )./range(ts,'all');
rescale_ts = @(ts, new_min, new_max) new_min + (new_max - new_min).*( ts - min(ts,[],'all') )./range(ts,'all');
