function dat = read_prog_input(inputfields, param_names, varargin)
% Originally for parsing input passed to ipcres program I created in Hitachi Cambridge 
% Laboratory, ~2010... So contains thing are not used in PtychoRunner, but does contain some useful
% functionality for checking function arguments. Should use Matlab inputParser functions in place of
% this though in the future.
%
% Arthur M. Blackburn, HCL, 2010 &
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

run_spec = varargin{1};

if ischar(run_spec)
    switch lower(run_spec)
        case 'direct'
            run_type = 1;
        case 'batch'
            run_type = 2;
        case 'struct'
            run_type = 3;
        case 'help'
            run_type = 4;
        otherwise
        error('Invalid run type specification. See help: "this_program.exe help."');
    end
else
    error('First argument must be run type specification. Too see help run "this_program.exe help".');
end

if run_type == 1
    if (length(varargin) >= 2) && ischar(varargin{2})  && ...
            isempty(str2num(varargin{2})) %#ok<ST2NM>
        % this last part checks that the input is not a char array
        % representing a number.
        dat.supplied_p.(inputfields(1,1)) = varargin{2};
    else
        error('Must supply a destination for the calculation output');
    end
    % read in the rest of the parameters supplied:
    if length(varargin) > 2
        for p_ind = 2:(length(varargin)-1)            
            field_name = inputfields{p_ind,1};                         
            dat.supplied_p.(field_name) = varargin{p_ind};                        
        end
    else
        error('Must supply at least one input');
    end
end


if run_type == 2
    % Collect cell array of left hand column names, searching for matches
    % from inputfields. When key phrase 'NEW RUN' is encountered, repeat
    % again and create a cell struct to contain the multiple runs. End the
    % reading with the phrase 'END RUN' is encountered.
    if (length(varargin) >= 2) && ischar(varargin{2}) && ...
            isempty(str2num(varargin{2})) %#ok<ST2NM>
        % this last part checks that the input is not a char array
        % representing a number.
        f_to_read = fopen(varargin{2});
    else
        error('Must supply a valid filename to read.');
    end
    file_cell = textscan(f_to_read,'%[^\n]');
    file_cell = file_cell{1};
    % get the end of run and end of file key word locations:
    run_markers = strcmpi(file_cell,'new run') | strcmpi(file_cell,'end run');
    run_ends = find(run_markers)-1;
    run_starts = [1; (run_ends(1:(end-1))+2) ];
    num_runs = length(run_ends);
    dat = cell(1,num_runs);
    for r_ind = 1:num_runs
        for l_ind = run_starts(r_ind):run_ends(r_ind)
            % look at each line for possible parameter matches:
            search_line = file_cell{l_ind};
            is_a_comment_line = sum(strcmp({'#','%'},search_line(1)));
            if ~is_a_comment_line
            for p_ind = 1:size(param_names,1)
                poss_names = param_names(p_ind,:);
                for poss_ind = 1:length(poss_names)
                    poss_entry = poss_names{poss_ind};
                    if ~isempty(poss_entry ) && ischar(poss_entry )
                        poss_name = poss_entry;
                    else
                        break;
                    end
                    find_ind = strfind(search_line, poss_name);
                    if ~isempty(find_ind)
                    % check that this text is either the last thing on the
                    % line (noting that the search line includes the CR),
                    % or has a space after it.                    
                    no_text_after = ...
                       ((find_ind+length(poss_name)) == length(search_line))  | ... 
                        strcmp( search_line(find_ind+length(poss_name)),' ');
                    % also check that there is nothing but whitespace
                    % before it:
                    no_text_before = ...
                        strcmp( search_line(1:(find_ind-1)), ...
                                char(double(' ')*ones(1,find_ind-1)));
                   
                    if  no_text_after && no_text_before
                                           
                        found_field = inputfields{p_ind,1}; 
                        % collect whatever is left over as a string:
                        dat{r_ind}.supplied_p.(found_field) = ...
                            search_line( (find_ind+length(poss_name)):(end-1));
                        % and there's no point looking any further
                        break;
                            
                    end
                    end
                end
            end
            end        
        end
    end    
end

if run_type == 3
    if (length(varargin) >= 2) && ischar(varargin{2})  && ...
            isempty(str2num(varargin{2})) %#ok<ST2NM>
        % then take this argument as a mat file structure containing the input structure:
        read_struct = load(varargin{2});
        read_name = fieldnames(read_struct); read_name = read_name{1};
        dat.supplied_p = read_struct.(read_name);
    else
        if isstruct(varargin{2})
            dat.supplied_p = varargin{2};
        else
            error('Must supply a structure as input');
        end
    end
end

if run_type == 4
    dat.help = true;
end

end

