% PLOT_STRUCTURE_FEATURE_DISTRIBUTIONS
% Adam Craig, 2023-06-26
% Plot the distributions of the 4 anatomical structural features.

hcp_data_header

subject_ids = training_subject_ids;
num_subjects = numel(subject_ids);
feature_matrices = cell(features_per_brain_area,1);
empty_feature_matrix = NaN(num_subjects, num_brain_areas);
for feature_index = 1:features_per_brain_area
    feature_matrices{feature_index} = empty_feature_matrix;
end
for subject_index = 1:num_subjects
    subject_id = subject_ids(subject_index);
    anat_file = get_brain_area_feature_file(subject_id);
    anat = load_brain_area_features(anat_file);
    for feature_index = 1:features_per_brain_area
        feature_matrices{feature_index}(subject_index,:) = anat(:,feature_index);
    end
end
for feature_index = 1:features_per_brain_area
    figure
    boxplot(feature_matrices{feature_index})
    xlabel('brain area')
    ylabel(brain_area_feature_names{feature_index})
end
for feature_index = 1:features_per_brain_area
    figure
    histogram( feature_matrices{feature_index}(:) )
    xlabel(brain_area_feature_names{feature_index})
    ylabel('count')
end
all_area_means = NaN(num_brain_areas, features_per_brain_area);
for feature_index = 1:features_per_brain_area
    area_mean = mean(feature_matrices{feature_index});
    all_area_means(:,feature_index) = area_mean;
    total_mean = mean( feature_matrices{feature_index}(:) );
    figure
    hold on
    histogram(area_mean)
    plot([total_mean total_mean],[0 num_brain_areas])
    legend({'area mean', 'total mean'})
    xlabel( sprintf('%s mean', brain_area_feature_names{feature_index}) )
    ylabel('count')
    ylim([0 num_brain_areas])
    hold off
end
write_binary_data(all_area_means, 'anatomy_mean.bin')
all_area_stds = NaN(num_brain_areas, features_per_brain_area);
for feature_index = 1:features_per_brain_area
    area_std = std(feature_matrices{feature_index});
    all_area_stds(:,feature_index) = area_std;
    total_std = std( feature_matrices{feature_index}(:) );
    figure
    hold on
    histogram(area_std)
    plot([total_std total_std],[0 num_brain_areas])
    legend({'area standard deviation', 'total standard deviation'})
    xlabel( sprintf('%s standard deviation', brain_area_feature_names{feature_index}) )
    ylabel('count')
    ylim([0 num_brain_areas])
    hold off
end
write_binary_data(all_area_stds, 'anatomy_std.bin')
for feature_index = 1:features_per_brain_area
    area_skewnesses = skewness(feature_matrices{feature_index});
    total_skewness = skewness( feature_matrices{feature_index}(:) );
    figure
    hold on
    histogram(area_skewnesses)
    plot([total_skewness total_skewness],[0 num_brain_areas])
    legend({'area skewnesses', 'total skewness'})
    xlabel( sprintf('%s kurtosis', brain_area_feature_names{feature_index}) )
    ylabel('count')
    ylim([0 num_brain_areas])
    hold off
end
for feature_index = 1:features_per_brain_area
    area_kurtosis = kurtosis(feature_matrices{feature_index});
    total_kurtosis = kurtosis( feature_matrices{feature_index}(:) );
    figure
    hold on
    histogram(area_kurtosis)
    plot([total_kurtosis total_kurtosis],[0 num_brain_areas])
    legend({'area kurtosis', 'total kurtosis'})
    xlabel( sprintf('%s kurtosis', brain_area_feature_names{feature_index}) )
    ylabel('count')
    ylim([0 num_brain_areas])
    hold off
end

