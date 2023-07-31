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

all_sc_mats = NaN(num_brain_areas, num_brain_areas, num_sc_training_subjects);
for s = 1:num_sc_training_subjects
    subject_id = sc_training_subject_ids(s);
    sc_file = get_structural_connectivity_file(subject_id);
    all_sc_mats(:,:,s) = load_structural_connectivity(sc_file);
end
sc_mean = mean(all_sc_mats,'all');
sc_std = std(all_sc_mats,0,'all');
fprintf('mean: %g, std: %g\n', sc_mean, sc_std)