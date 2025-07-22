function probe_positions_0 = gen_coords_for_ptychoshelves(source_file_or_data, output_name, nm_V, rot_angle_deg, convert_to_angstrom, ....
    do_transpose)
% A function to create a position file in a format PtychoShelves prefers from an Azorus hpl file, or
% similar. The main is that I have 'coded' the transformation case numbers here (1 - 8) to deal with
% all the combinations of rotate, mirror, flip etc. Dealing explicitly with rotate, flip, mirror can
% be confusing, as it depends on the order of operations, so I went for transformation number which
% removes the worry about the ordering of rotates, flips, transposes etc.
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

if nargin < 6
   do_transpose = 3; % this does nothing x,y -> x,y
end

% do_transpose = 4; % this works for the data I acquired on 2020/Dec/24.
% Here note the SU9000 view was transposed (x and y)... In 'early days' check Azorus config files, 
% and that no-one has messed with x, y cables between runs or when you weren't looking!

if nargin < 5
    convert_to_angstrom = true;
end

if convert_to_angstrom
    A_nm = 10;
else
    A_nm = 1;
end


if ischar(source_file_or_data)
    [basedir, ff, ext] = fileparts(source_file_or_data);
    fname = [ff, ext];
    if isempty(output_name)
        outname = [ff,'_CONV.hdf5'];
        output_name = fullfile(basedir,outname);
    end
    loadup = fullfile(basedir, fname);
    exp_data = loadHPL(loadup);
    crds.x_arr = exp_data.coords(:,1)';
    crds.y_arr = exp_data.coords(:,2)';
else
%     if isempty(output_name)
%         error('Output name could not be determined from input');        
%     end
    src_shape = size(source_file_or_data);
    if length(src_shape) == 2
        small_ax = find(src_shape == 2);
        if length(small_ax) > 1
            error('Cannot determine x,y axes of input coord data.');
        else
            if small_ax == 1
                crds.x_arr = source_file_or_data(1,:);
                crds.y_arr = source_file_or_data(2,:);
                
            elseif small_ax == 2
                crds.x_arr = source_file_or_data(:,1)';
                crds.y_arr = source_file_or_data(:,2)';
                
            end
        end
    else
        error('Shape of input data is incorrect.');
    end       
end

%%
% Just check transpose was used by comparing plots here:
% figure; 
% tiledlayout(2,4)
% for trans_coords = 1:8
    % Apparently number 4 does the trick...
% for trans_coords = 4

        switch(do_transpose)
            case 1
                exp_data.coords = [-crds.y_arr',  crds.x_arr'];
                fprintf('Transposing coordinates case: 1!\n');
            case 2
                exp_data.coords = [ crds.y_arr',  crds.x_arr']; %
                fprintf('Transposing coordinates case: 2!\n');
            case 3
                exp_data.coords = [ crds.x_arr',  crds.y_arr']; %
                fprintf('Transposing coordinates case: 3!\n');
            case 4 
                exp_data.coords = [ crds.x_arr', -crds.y_arr']; % 
                fprintf('Transposing coordinates case: 4!\n');
            case 5 
                exp_data.coords = [ crds.y_arr', -crds.x_arr']; % 
                fprintf('Transposing coordinates case: 5!\n');
            case 6
                exp_data.coords = [-crds.y_arr', -crds.x_arr']; % 
                fprintf('Transposing coordinates case: 6!\n');
            case 7 
                exp_data.coords = [-crds.x_arr', crds.y_arr']; % 
                fprintf('Transposing coordinates case: 7!\n');
            case 8 
                exp_data.coords = [-crds.x_arr', -crds.y_arr']; % 
                fprintf('Transposing coordinates case: 8!\n');
            otherwise
                error('Invalid case # for transforming coords');
        end
        % Do the rotation here:
        rot_angle = (rot_angle_deg /180)*pi;  

        R = [cos(rot_angle) -sin(rot_angle);...
             sin(rot_angle)  cos(rot_angle)];
            
        xy_posns = exp_data.coords.';
        rotated_posns = R*xy_posns;
        exp_data.coords = rotated_posns.';
        %
        

%%
probe_positions_0 = exp_data.coords*nm_V*A_nm;
% data_file_outstring = 'data_postions_%02.0f';        
% savename = sprintf(data_file_outstring, trans_coords);
if ~isempty(output_name)
    [~, ~, crds_EXT] = fileparts(output_name);
    if strcmpi(crds_EXT,'.hdf5')
        hdf5write(output_name,'/probe_positions_0',probe_positions_0);
    elseif strcmpi(crds_EXT,'.mat')
        save(output_name,'probe_positions_0','-v7.3');
    else
        error('Unregonized Poistion File Output Format');
    end

    fprintf('Positions saved to file: %s\n',output_name);
end

%% Test read:
%{
figure;
pos_in_file = loadHPL(fullsavename);
figure; plot(pos_in_file.probe_positions_0(1,:),pos_in_file.probe_positions_0(2,:)); title(sprintf('As in %s',savename),'interpreter','none');
axis image
hold on; plot(pos_in_file.probe_positions_0(1,1:10),pos_in_file.probe_positions_0(2,1:10),'ro');
%}
end



