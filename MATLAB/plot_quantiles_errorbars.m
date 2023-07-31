function plot_quantiles_errorbars(values_table, groupvar, datavar, low_quant, mid_quant, high_quant, linspec, max_num_points)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if ~exist('groupvar','var')
    groupvar = 'epoch';
end
if ~exist('datavar','var')
    datavar = 'loss';
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
if ~exist('linspec','var')
    linspec = '-r';
end
if ~exist('max_num_points','var')
    max_num_points = 30;
end
quantiles_table = grpstats(values_table, groupvar, {@(x) quantile(x,low_quant) @(x) quantile(x,mid_quant) @(x) quantile(x,high_quant)}, 'DataVars', {datavar}, 'VarNames', {'GrpVar', 'GrpCount', 'LowQuant', 'MidQuant', 'HighQuant'});
num_points = size(quantiles_table,1);
if num_points > max_num_points
    selected_indices = floor( linspace(1,num_points,max_num_points) );
    quantiles_table = quantiles_table(selected_indices,:);
end
errorbar(quantiles_table.GrpVar, quantiles_table.MidQuant, quantiles_table.MidQuant - quantiles_table.LowQuant, quantiles_table.HighQuant - quantiles_table.MidQuant, linspec)

end