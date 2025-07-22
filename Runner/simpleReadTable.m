function  table_data = simpleReadTable( fname, sheet_name, region)
% A simple function for reading in excel tables
%    simpleReadTable( fname, sheet_name, region)
%
% fname - file name (*.xls, etc);
% sheetname - excel sheet name
% region - [first row number, first column number, number of rows, number
% of columns]
% 
% function will expect:
%  - Headers to be in first row,
%  - Data types to be in second row, e.g. double, string, char, etc.
%  - Data to follow in the third row and or below.
%
% The function will read  columns upto max of 'number of columns', BUT as
% sson as column is found that contains no header information (i.e. is empty),
% the reading stops. Thus, make sure there are no blank headers between data columns!
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

