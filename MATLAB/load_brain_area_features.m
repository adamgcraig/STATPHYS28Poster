function brain_area_features = load_brain_area_features(brain_area_feature_file, num_brain_areas, features_per_brain_area)
%LOAD_BRAIN_AREA_FEATURES Load HCP structual MRI area features from a binary file.
%   We need to specify the dimensions,
%   because we store the data as an array of 64-bit floats in a binary.

if ~exist('num_brain_areas','var')
    num_brain_areas = 360;
end
if ~exist('features_per_brain_area','var')
    features_per_brain_area = 4;
end
brain_area_features = load_data_from_binary(brain_area_feature_file, num_brain_areas, features_per_brain_area);

end