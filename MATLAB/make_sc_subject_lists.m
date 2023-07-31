hcp_data_header

name = {'training' 'validation' 'testing'};
list = {training_subject_ids, validation_subject_ids, testing_subject_ids};
num_lists = numel(list);
for list_index = 1:num_lists
    list_name = name{list_index};
    subject_ids = list{list_index};
    num_ids = numel(subject_ids);
    has_sc_data = false(num_ids,1);
    for s = 1:num_ids
        subject_id = subject_ids(s);
        sc_file = get_structural_connectivity_file(subject_id);
        has_sc_data(s) = exist(sc_file,'file') > 0;
    end
    sc_subject_ids = subject_ids(has_sc_data);
    sc_subject_id_strings = compose('%u', sc_subject_ids);
    writelines( sc_subject_id_strings, [data_dir filesep sprintf('sc_%s_subject_ids.txt',list_name)] )
end