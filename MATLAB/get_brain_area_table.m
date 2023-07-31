function brain_area_table = get_brain_area_table(brain_area_table_file)
%GET_BRAIN_AREA_TABLE Load table of brain area information.
%   brain_area_table_file: See MMP360coordinator.xlsx.
%   brain_area_table: has 4 columns
%   x, y, and z are the coordinates of the centroid of each area.
%   name is the area name (abreviation).
%   The L or R at the start of the name indicates left or right hemisphere.

brain_area_cell = readcell(brain_area_table_file);
% The first row is just the word "#mmp360".
x = brain_area_cell(2:end,1);
y = brain_area_cell(2:end,2);
z = brain_area_cell(2:end,3);
% For some reason,
% columns 4 and 5 are just the numbers 1 & 2, respectively.
name = brain_area_cell(2:end,6);
brain_area_table = table(x,y,z,name);

end