function [experi_def, set_params, recipe_variations] = LoadFromSpreadSheet(FileLocation, ...
    ExcelFileName, ExperimentSheet, RecipeSheet, RunOnSamples)
% Loads data for ptychography reconstructions from spreadsheets
%
% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2021
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
% 
% *******************************************************************

first_col = 2; num_cols = 100; 
end_col = first_col + num_cols  - 1;
% convert to Excels nasty format:
first_col_A1 = idx2A1(first_col); 
end_col_A1   = idx2A1(end_col);

first_row = 2; num_rows = 1000; 
end_row = first_row + num_rows - 1;

experi_def.localsummary_dir = FileLocation;
experi_def.summary_table.xl_filename = ExcelFileName;
% if exist('log_file','var')
%     log_filename = fullfile(experi_def.localsummary_dir, log_file);
%     diary(log_filename);
%     fprintf('Execution Started at %s\n',datestr(datetime('now')));
% end
experi_def.summary_table.sheet_name = ExperimentSheet;
experi_def.summary_table.xlsheet_range = sprintf('%s%.0f:%s%.0f',first_col_A1,first_row, end_col_A1,end_row);
experi_def.summary_table.run_on_samples = RunOnSamples;

experi_def.summary_table.use_ys_to_define_experis = true;
experi_def.summary_table.load_from_mat = true;
experi_def.summary_table.use_spreadsheet_experi_params = true;
experi_def.summary_table.load_from_mat_when_possible = true;% Use a subregion?

experi_def.recipe_table.xl_filename = ExcelFileName;
% experi_def.recipe_table.sheet_name = 'Recipes'; %as used for ePIE on
% Phasefocus box
experi_def.recipe_table.sheet_name = RecipeSheet;
experi_def.recipe_table.xlsheet_range = sprintf('%s%.0f:%s%.0f',first_col_A1,first_row,end_col_A1,end_row);
experi_def.recipe_table.prime_probe_num = 1;

% collect the header information, which includes the types of the data!
% summary_opts_hdr = spreadsheetImportOptions("NumVariables", num_cols);
% summary_opts_hdr.Sheet = experi_def.summary_table.sheet_name;
% summary_opts_hdr.DataRange = sprintf('%s%.0f:%s%.0f',first_col_A1,first_row,end_col_A1,first_row+1);
% summary_opts_hdr.ReadVariableNames = true;

%{
set_param_headers = readtable(fullfile(experi_def.localsummary_dir, experi_def.summary_table.xl_filename), ...
    summary_opts_hdr, 'ReadVariableNames', true);
% the data types are in row 2 of the above. For some reason the above also
% ignores 'ReadVariableNames', true

% find the actual number of used columns from the above;
ii = 1;
while ii <= num_cols
    if isempty(set_param_headers{1,ii}{:})
        break
    end
    ii = ii + 1;
end
used_cols = ii - 1;

set_params_opts = spreadsheetImportOptions("NumVariables", used_cols);
set_params_opts.VariableNames = set_param_headers{1,1:used_cols};
set_params_opts.VariableTypes = set_param_headers{2,1:used_cols};
set_params_opts.DataRange = sprintf('%s%.0f:%s%.0f', first_col_A1, first_row+2, idx2A1(first_col + used_cols - 1), end_row);
set_params_opts.Sheet = experi_def.summary_table.sheet_name;
%}

fname = fullfile(experi_def.localsummary_dir, experi_def.summary_table.xl_filename);

set_params = myreadtable(fname, ...
    experi_def.summary_table.sheet_name,...
    [first_row, first_col, num_rows, num_cols]);

% set_params = readtable(fullfile(experi_def.localsummary_dir, experi_def.summary_table.xl_filename), ...
%     set_params_opts);

% set_params = readtable(fullfile(experi_def.localsummary_dir, experi_def.summary_table.xl_filename), ...
%     set_params_opts, ...
%     'Sheet', experi_def.summary_table.sheet_name,'Range', experi_def.summary_table.xlsheet_range, 'ReadVariableNames',true);


recipe_variations = myreadtable( fname, ...
    experi_def.recipe_table.sheet_name, ...
    [first_row, first_col, num_rows, num_cols]);


% Read in info on the recipe variations to be tried;
% recipe_variations = readtable(fullfile(experi_def.localsummary_dir, experi_def.recipe_table.xl_filename), ...
%     'Sheet', experi_def.recipe_table.sheet_name,'Range', experi_def.recipe_table.xlsheet_range, 'ReadVariableNames',true);

end


function table_data = myreadtable( fname, sheet_name, region)

first_row = region(1);
first_col = region(2);
first_col_A1 = idx2A1(first_col);

num_rows = region(3);
num_cols = region(4);

end_col_A1 = idx2A1(first_col + num_cols - 1);
end_row = first_row + num_rows - 1;

% Collect the header information, which includes the types of the data!
hdr_opts = spreadsheetImportOptions("NumVariables", num_cols);
hdr_opts.Sheet = sheet_name;
hdr_opts.DataRange = sprintf('%s%.0f:%s%.0f', first_col_A1, first_row, end_col_A1, first_row+1);

headers = readtable(fname, hdr_opts, 'ReadVariableNames', true);

% Find the actual number of used columns from the above;
ii = 1;
while ii <= num_cols
    if isempty(headers{1,ii}{:})
        break
    end
    ii = ii + 1;
end
used_cols = ii - 1;

table_opts = spreadsheetImportOptions("NumVariables", used_cols);
table_opts.VariableNames = headers{1,1:used_cols};
table_opts.VariableTypes = headers{2,1:used_cols};
table_opts.DataRange = sprintf('%s%.0f:%s%.0f', first_col_A1, first_row+2, idx2A1(first_col + used_cols - 1), end_row);
table_opts.Sheet = sheet_name;

table_data = readtable(fname, table_opts);

end
