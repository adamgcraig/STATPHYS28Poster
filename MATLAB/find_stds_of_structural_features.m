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

all_anat_mats = NaN(num_brain_areas, features_per_brain_area, num_sc_training_subjects);
for s = 1:num_sc_training_subjects
    subject_id = sc_training_subject_ids(s);
    anat_file = get_brain_area_feature_file(subject_id);
    all_anat_mats(:,:,s) = load_brain_area_features(anat_file);
end
% sc_mask = sum(all_anat_mats > 0, 3)/num_sc_training_subjects;
% num_non0_pairs = nnz(sc_mask);
% num_pairs = numel(sc_mask);
% write_binary_data(sc_mask,'E:\HCP_data\dtMRI_binaries\sc_mask.bin');
% fprintf('%u of %u area pairs are non-0 at least once.\n', num_non0_pairs, num_pairs);
% figure
% histogram( sc_mask(:) )
% xlabel('fraction of training subjects with a non-0 SC')
% ylabel('number of brain area pairs')
thickness = all_anat_mats(:,1,:);
myelination = all_anat_mats(:,2,:);
curvature = all_anat_mats(:,3,:);
sulcus_depth = all_anat_mats(:,4,:);

disp('means')
thickness_mean = mean(thickness,'all');
myelination_mean = mean(myelination,'all');
curvature_mean = mean(curvature,'all');
sulcus_depth_mean = mean(sulcus_depth,'all');
all_mean = mean(all_anat_mats,'all');
fprintf('thickness: %f\n',thickness_mean)
fprintf('myelination: %f\n', myelination_mean)
fprintf('curvature: %f\n', curvature_mean)
fprintf('sulcus depth: %f\n', sulcus_depth_mean)
fprintf('all structural features: %f\n', all_mean)

disp('stds')
thickness_std = std(thickness,0,'all');
myelination_std = std(myelination,0,'all');
curvature_std = std(curvature,0,'all');
sulcus_depth_std = std(sulcus_depth,0,'all');
all_std = std(all_anat_mats,0,'all');
fprintf('thickness: %f\n',thickness_std)
fprintf('myelination: %f\n', myelination_std)
fprintf('curvature: %f\n', curvature_std)
fprintf('sulcus depth: %f\n', sulcus_depth_std)
fprintf('all structural features: %f\n', all_std)
