function create_ptyshelves_subset(dp_file_in, pos_file_in, ...
                                  dp_file_out, pos_file_out, ...
                                  shape_in, subregion_out,...
                                  nm_V, rot_angle_deg, do_transpose, ...
                                  mask_file, append_probe)
% Function to create a region subset for pytchoshelves
% Can also use this to do a straigh conversion if subregion_out = [];
%
% create_ptyshelves_subset(dp_file_in, pos_file_in, ...
%                          dp_file_out, pos_file_out, ...
%                          shape_in, subregion_out,...
%                          nm_V, rot_ang_deg)
% 
% dp_file_in:  can be mat,hp, or hpl file.
% pos_file_in: need not be specified if positions are in the 'dp_file', but
% note if using hpl file, then nm_V and rot_ang_deg must be specified to
% create file to be used for ptychoshevles_set.
% shape_in: is the shape of the mesh used to define the DP acquire
% positions. This can be left empty if this data is in the HPL file in
% exp_data.attrs.meshParams.shape, where exp_data is the root of the HPL file;
% 
% If xxx_file_out is left empty default names of 
%   ...SUBSET_dp_data.mat and 
%   ...SUBSET_data_position.hdf5
% will be used
%
% For example: 
% 
% create_ptyshelves_subset('D:\Data\Arthur\20201224\data_dp.mat',...
%                          'D:\Data\Arthur\20201224\data_position.hdf5',...
%                          'D:\Data\Arthur\20201224\data_dp_30by30.mat',...
%                          'D:\Data\Arthur\20201224\data_position_30by30.hdf5',...
%                          [45, 45], [30, 30]);
%
% Also used with "...\Runner\TestScripts\Part_0_Prepare_Data.m"
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


if nargin< 7
    nm_V = 1;
end

if nargin < 8
    rot_angle_deg = 0;
end

if nargin < 9
    do_transpose = 3;
end

if nargin >= 10
    if ~isempty(mask_file)
        mask = load(mask_file);
        hasmask = true;
    else
        hasmask = false;
    end
else
    hasmask = false;
end

if nargin < 11
    append_probe = false;
end


[basedir, fname, ext] = fileparts(dp_file_in);
if strcmpi(ext,'.mat')
    % Note: section below is also in L:\em_utils\Ptycho\Processing\PtychoProcessFromTableF_MS.m
    dp_set = load(dp_file_in);
   
    if isfield(dp_set,'dat') % then has come from ?
%           if isfield(dp_set,'diffPats') % then has come from oops?.
        permute_order_of_loaded_dps = [3 2 1];
        dp = permute(dp_set.dat, permute_order_of_loaded_dps);
    elseif isfield(dp_set,'dps') % then this has come from my multislice simulation, and is same order as from loadHPL
        permute_order_of_loaded_dps = [2 3 1]; % ExamplePrepareData saves as name 'dps', so this is taken as order. Same as from loadHPL...
        dp = permute(dp_set.dps, permute_order_of_loaded_dps);
    elseif isfield(dp_set, 'normed_images')
        dp = zeros([size(dp_set.normed_images{1}),length(dp_set.normed_images)]);
        permute_order_of_loaded_dps = [1 2 3]; % no need to permute in this case.
        for ii=1:length(dp_set.normed_images)
            dp(:,:,ii) = dp_set.normed_images{ii};
        end
    elseif isfield(dp_set, 'dp') % then has come from Chen Ptychoset EM, or this file!
        permute_order_of_loaded_dps = [1 2 3];
        dp = dp_set.dp;
        % for now this is the only one that supports embedded mask

    else
        error('Unregonized mat file data structure.');
    end
    if isfield(dp_set, 'probe_0')
        probe0 = dp_set.probe_0;
    end
elseif strcmpi(ext,'.hpl') || strcmpi(ext,'.hp')
    exp_data = loadHPL(dp_file_in);
    if ~isfield(exp_data,'diffPats')
        if isfield(exp_data,'diffraction')
            if isfield(exp_data.diffraction,'micrograph')
                dp = exp_data.diffraction.micrograph; % when loaded has something like:     micrograph: [2025×512×512 single]
                exp_data.diffraction = rmfield(exp_data.diffraction,'micrograph');
                permute_order_of_loaded_dps = [2, 3, 1];
                dp = permute(dp, permute_order_of_loaded_dps);    % this makes it look like as seen in azorus, but with index as last element.
            else
                error('Cannot find diffraction data in loaded HP file');
            end                       
        end
        if (isfield(exp_data,'attrs') && isfield(exp_data.attrs,'meshParams') && isfield(exp_data.attrs.meshParams,'shape'))
            shp_in_file = exp_data.attrs.meshParams.shape + 1; % currently Rob stores this as size - 1
            if ~isempty(shape_in) && any(shape_in ~= shp_in_file)
                warning('Shape in file is used in precedence to that supplied as function input.');
                shape_in = shp_in_file;
            end
        end
    else
        dp = exp_data.diffPats; % for Robb's older hpl files...
%         permute_order_of_loaded_dps = [2, 1, 3]; % might need to swap first 2 elements...
        permute_order_of_loaded_dps = [1, 2, 3]; % might need to swap first 2 elements...
        dp = permute(dp, permute_order_of_loaded_dps); 
    end
elseif isempty(ext) && strcmpi(fname,'none')
    dp = [];
    hasmask = false;
else
    error('Unknown source file for diffraction data.');
end

if hasmask
    warning('Assuming that loaded mask has same ordering as loaded dp_data');
    mask = permute(mask.mask, permute_order_of_loaded_dps); % quick hack.... FIX - TO DO
end



if isempty(pos_file_in) 
    if exist('exp_data','var') % Big nasty shortcut, here if exp_data is created, then a hpl must have be been read (we assume), and it shoudl contiain 'coords' field....
        crds.x_arr = exp_data.coords(:,1)';
        crds.y_arr = exp_data.coords(:,2)'; % would need to to do conversion to nm, A etc, if we're actually going to apply this...
        convert_to_angstrom = true;
    else
%         % Are the coords within the HPL file?
%         if strcmpi(ext,'.hpl') || strcmpi( ext, '.hp')
%             
%             convert_to_angstrom = true;
%         else
%         end
    end
    % it appears that when coming from the old hpl, the coords are
    % (mostly...) transposed.
%     if strcmpi(ext,'.hpl')           
%         crds.x_arr = exp_data.coords(:,2)';
%         crds.y_arr = exp_data.coords(:,1)'; % 
%     end
    
else
    [~, ~, crds_EXT] = fileparts(pos_file_in);
    if strcmpi(crds_EXT,'.hdf5')
        posns = h5read(pos_file_in,'/probe_positions_0');            
        crds.x_arr  = posns(:,1).'; 
        crds.y_arr  = posns(:,2).'; 
        if exist('exp_data','var') && isfield(exp_data,'coords')
            warning('Coords data is supplied in source DP hp/hpl file, but position data from independent hdf5 file was used instead.');
        end
        convert_to_angstrom = false; % as data is already in angstrom if read from this format.
    end
    if strcmpi(crds_EXT,'.mat')
        if ~strcmpi(pos_file_in, dp_file_in) 
            dp_set = load(pos_file_in);       
        end
        if isfield(dp_set,'voltage_coords')
            warning('Carefully check voltage_coords section out as it has not been tested.\n');
            crds.x_arr = dp_set.voltage_coords(:,1)'*nm_V;
            crds.y_arr = dp_set.voltage_coords(:,2)'*nm_V;                        
        elseif isfield(dp_set,'probe_0') % then this has come from my multislice simulation, where 1 unit is 1A, ie "nm_per_V is 0.1"
%             exp_data.coords = [crds.y_arr',  -crds.x_arr']; % These are actual positions. #5 - the swap later will make everything as normal.
% 1:(+x,+y); 2:+x,-y; 3:-x,+y; 4:-x,-y; 5:-y,-x; 6:-y,x; 7: y,-x; 8:y, x
            crds.x_arr =  dp_set.y_arr;
            crds.y_arr =  dp_set.x_arr;
            % also load the probe model at this point:
%             exp_data.probeModel = crds.probe_0;
        elseif isfield(dp_set,'req_postions')
            crds.x_arr = crds.req_postions(1:(end-1),1)';
            crds.y_arr = crds.req_postions(1:(end-1),2)';
%             exp_data.coords = [crds.x_arr' , crds.y_arr'];
        else
            error('Unregonized mat file data structure.');
        end
        convert_to_angstrom = true;
    end

end  
    
ps = ptychoset_mult;
coords_name = 'read_in_coords';
ps.dps.attrs.meshParams.shape = shape_in;

ps.dps.( coords_name)(1,:) = crds.x_arr; % xpts
ps.dps.( coords_name)(2,:) = crds.y_arr; % ypts
 
if ~isempty(subregion_out) && ~(all(subregion_out == shape_in)) 
    [dp_posns, dp_indices] = ps.create_subset(subregion_out, coords_name);
    dp = dp(:,:,dp_indices);
    if hasmask
        mask = mask(:,:,dp_indices);
    end
else
    dp_posns = [crds.x_arr; crds.y_arr];
end

if isempty(dp_file_out)
    if ~hasmask
        dp_file_out = fullfile(basedir, [fname,'_SUBSET_data_dp.mat']);
    else
        dp_file_out = fullfile(basedir, [fname,'_SUBSET_mask_data_dp.mat']);
    end
end

if ~isempty(dp)
    if ~hasmask
        save(dp_file_out, 'dp','-v7.3');
        fprintf('DPs saved into file %s\n',dp_file_out);
    else
        save(dp_file_out, 'dp','mask','-v7.3');
        fprintf('DPs saved into file with mask embedded %s\n',dp_file_out);
    end
end

if isempty(pos_file_out)
    pos_file_out = fullfile(basedir, [fname,'_SUBSET_data_position.hdf5']);
end

gen_coords_for_ptychoshelves(dp_posns, pos_file_out, nm_V, rot_angle_deg, convert_to_angstrom, do_transpose);

if ~isempty(dp)
    if append_probe
        if exist('probe0','var')
            save(dp_file_out,'probe0','-append');
        end
    end
end
end