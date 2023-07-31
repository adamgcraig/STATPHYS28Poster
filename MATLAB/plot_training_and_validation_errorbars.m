function fig = plot_training_and_validation_errorbars(training_table, validation_table, loss_var, epoch_var, low_quant, mid_quant, high_quant)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if ~exist('epoch_var','var')
    epoch_var = 'epoch';
end
if ~exist('loss_var','var')
    loss_var = 'loss';
end
if ~exist('low_quant','var')
    low_quant = 0.05;
end
if ~exist('mid_quant','var')
    mid_quant = 0.50;
end
if ~exist('high_quant','var')
    high_quant = 0.95;
end
fig = figure;
hold on
plot_quantiles_errorbars(training_table, epoch_var, loss_var, low_quant, mid_quant, high_quant, '-r')
plot_quantiles_errorbars(validation_table, epoch_var, loss_var, low_quant, mid_quant, high_quant, '--g')
legend({'training', 'validation'})
xlabel(epoch_var)
ylabel(  sprintf( '%.0f-th, %.0f-th, %.0f-th percentiles of %s', 100*low_quant, 100*mid_quant, 100*high_quant, strrep(loss_var,'_',' ') )  )
hold off
end