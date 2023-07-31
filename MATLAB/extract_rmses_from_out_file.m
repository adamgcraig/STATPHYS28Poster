% out_file_path = 'C:\Users\agcraig\Documents\GitHub\IzhikevichReservoirComputer\segmented_reservoir_computing_leave_one_out_test_23089.out';
out_file_path = 'C:\Users\agcraig\Documents\GitHub\IzhikevichReservoirComputer\segmented_reservoir_computing_leave_one_out_test_23089.out';
out_file_str = readlines(out_file_path);

training_rep_expression = 'time (?<time>\d+\.\d) seconds, time series (?<subject>\d+) (?<suffix>\d_\w\w), training rep (?<rep>\d+), RMSE (?<tsrmse>\d\.\d+), FC RMSE (?<fcrmse>\d\.\d+)';
training_tokens = regexp(out_file_str,training_rep_expression,'names');
training_rmse_table = struct2table([training_tokens{:}]);
training_rmse_table.time = str2double(training_rmse_table.time);
training_rmse_table.subject = str2double(training_rmse_table.subject);
training_rmse_table.rep = str2double(training_rmse_table.rep);
training_rmse_table.tsrmse = str2double(training_rmse_table.tsrmse);
training_rmse_table.fcrmse = str2double(training_rmse_table.fcrmse);
disp('training results')
disp( training_rmse_table(1:10,:) )

validation_rep_expression = 'time (?<time>\d+\.\d) seconds, time series (?<subject>\d+) (?<suffix>\d_\w\w), validation rep, RMSE (?<tsrmse>\d\.\d+), FC RMSE (?<fcrmse>\d\.\d+)';
validation_tokens = regexp(out_file_str,validation_rep_expression,'names');
validation_rmse_table = struct2table([validation_tokens{:}]);
validation_rmse_table.time = str2double(validation_rmse_table.time);
validation_rmse_table.subject = str2double(validation_rmse_table.subject);
validation_rmse_table.tsrmse = str2double(validation_rmse_table.tsrmse);
validation_rmse_table.fcrmse = str2double(validation_rmse_table.fcrmse);
disp('validation results')
disp( validation_rmse_table(1:10,:) )

combined_ts_rmses = [ training_rmse_table.tsrmse; validation_rmse_table.tsrmse ];
combined_fc_rmses = [ training_rmse_table.fcrmse; validation_rmse_table.fcrmse ];
combined_labels = [
    repmat( {'training'}, size(training_rmse_table,1), 1 );
    repmat( {'validation'}, size(validation_rmse_table,1), 1 )
    ];
figure
boxplot(combined_ts_rmses, combined_labels)
ylabel('time series RMSE')
figure
boxplot(combined_fc_rmses, combined_labels)
ylabel('FC RMSE')
