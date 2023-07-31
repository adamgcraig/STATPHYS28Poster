function structural_connectivity = load_structural_connectivity(structural_connectivity_file, num_brain_areas)
%LOAD_STRUCTURAL_CONNECTIVITY Load HCP DT-MRI structural connectivity matrix from a binary file.
%   We need to specify the dimensions,
%   because we store the data as an array of 64-bit floats in a binary.

if ~exist('num_brain_areas','var')
    num_brain_areas = 360;
end
structural_connectivity = load_data_from_binary(structural_connectivity_file, num_brain_areas, num_brain_areas);

end