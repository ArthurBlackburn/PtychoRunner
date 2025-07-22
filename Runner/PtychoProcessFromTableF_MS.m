function PS_recon = PtychoProcessFromTableF_MS(exec_params, experi_def, set_params, run_inds, recipe_variations)
% This function is for processing reconstructions guided by parameters in
% Excel sheets. This "_MS" version works with the multislice variant of ptychoshelves
% as developed by Z Chen et al.
%
% It is called from RunReconN. Examples of RunRecon are given in 
% ...\TestScripts\Part_1and2_Reconstructions.m
%
%% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2020
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
% 
% *******************************************************************

try

start_time = now;
%%

%% Loop setup
% ext = [];
testing = false;
% exec_params.use_pi_box = false;
% just a couple of handy helper functions:
contains_num = @(x) all(~isempty(x) & isnumeric(x) & ~isnan(x));
cell_is_true = @(x) ~isempty(x) && ~isa(x,'string') && ~isnan(x) && logical(x); 
% cell_is_char = @(x) (~isempty(x) && ~any(isnan(x)) && ischar(x));
cell_is_charorstring = @(x) ~isempty(x) && (isa(x,'string')||isa(x,'char'));

% if exec_params.show_last_recon
%     f1 = figure;
%     last_montage = axes;
%     set(f1,'Name', 'Last Recon');
% end

%%
% for set_number = run_inds(1)
for set_number = run_inds
    
    subregion = [];
    if any(strcmpi(fieldnames(experi_def), 'subregion_override'))
        subregion = experi_def.subregion_override;
    else
        if any(strcmpi(set_params.Properties.VariableNames, 'subregion'))
            if ~isnan(set_params.subregion(set_number))
                subregion = logical(set_params.subregion(set_number));
            else
                subregion = [];
            end
        else
            subregion = [];
        end
    end
    
%%    
    base_dir = set_params.(exec_params.base_dir_fieldname){set_number};
    if any(strcmpi(set_params.Properties.VariableNames, 'recon_savepath'))
        reconsavepath = set_params.recon_savepath{set_number};
    else
        reconsavepath = base_dir;
    end
    
    fname = set_params.file_name{set_number};
    if any(strcmpi(set_params.Properties.VariableNames, 'extension'))
        ext = set_params.extension{set_number};
        fname_base = fname;
    else
        [~,fname_base,ext] = fileparts(fname);
        ext = ext(2:end); % remove the dot.
    end
    loadup = fullfile(base_dir, [fname_base,'.', ext]);
    if any(strcmpi(set_params.Properties.VariableNames, 'trans_coords'))
        if ~isnan(set_params.trans_coords(set_number))
            trans_coords = set_params.trans_coords(set_number);
        else
            % should we set a default transformation??
            set_params.trans_coords(set_number) = 3;
            trans_coords = set_params.trans_coords(set_number); % This leaves it IN:[X,Y] -> OUT:[X, Y]
        end
    end
    save_name = set_params.save_name{set_number};
    if exist('crds','var') % should really sort shit out to remove this nasty global.
        clear crds
    end
    
    if experi_def.summary_table.load_from_mat_when_possible && strcmpi(ext,'mat')
        experi_def.summary_table.load_from_mat = true;
    else
        experi_def.summary_table.load_from_mat = false;
    end
    % check that incompatible parameters have not been given.
    if experi_def.summary_table.load_from_mat 
        experi_def.summary_table.use_spreadsheet_experi_params = true; 
    end
    
    if ~testing
        if experi_def.summary_table.load_from_mat
            dp_set = load(loadup);

            if isfield(dp_set,'dat') % then has come from ?
%           if isfield(dp_set,'diffPats') % then has come from oops?.
                permute_order_of_loaded_dps = [3 2 1];
                exp_data.diffPats = permute(dp_set.dat, permute_order_of_loaded_dps);
            elseif isfield(dp_set,'dps') % then this has come from my multislice simulation, and is same order as from loadHPL
                permute_order_of_loaded_dps = [2 3 1];
                exp_data.diffPats = permute(dp_set.dps, permute_order_of_loaded_dps);
            elseif isfield(dp_set, 'normed_images')
                exp_data.diffPats = zeros([size(dp_set.normed_images{1}),length(dp_set.normed_images)]);
                permute_order_of_loaded_dps = []; % no need to permute in this case.
                for ii=1:length(dp_set.normed_images)
                    exp_data.diffPats(:,:,ii) = dp_set.normed_images{ii};
                end
            elseif isfield(dp_set, 'dp') % then has come from Chen Ptychoset EM
                permute_order_of_loaded_dps = [1 2 3];
                exp_data.diffPats = dp_set.dp;                
            else
                error('Unregonized mat file data structure.');
            end
            
            if isfield(dp_set,'probe0')
                exp_data.probe0 = dp_set.probe0;
            end
            
            if any(strcmpi(set_params.Properties.VariableNames, 'coords'))
                [crds_PATHSTR,crds_NAME,crds_EXT] = fileparts(set_params.coords{set_number});
                if isempty(crds_PATHSTR) || strcmpi('..',crds_PATHSTR) ||  strcmpi('.',crds_PATHSTR)
                    crds_file_name = fullfile(base_dir, set_params.coords{set_number});
                else
                    crds_file_name = set_params.coords{set_number};
                end
                
                if strcmpi(crds_EXT,'.mat')
                    if ~strcmpi(crds_file_name, loadup) 
                        crds = load(crds_file_name);
                    else
                        crds = dp_set;
                    end
                    if isfield(crds,'voltage_coords')
                        exp_data.coords = crds.voltage_coords;
                        crds.x_arr = exp_data.coords(:,1)';
                        crds.y_arr = exp_data.coords(:,2)';
                        collection_params = [];
                    elseif isfield(crds,'probe_0') % then this has come from my multislice simulation, where 1 unit is 1A, ie "nm_per_V is 0.1"
                        exp_data.coords = [crds.y_arr',  -crds.x_arr']; % These are actual positions. #5 - the swap later will make everything as normal.

                        % also load the probe model at this point:
                        exp_data.probeModel = crds.probe_0;
                    elseif isfield(crds,'req_postions')
                        crds.x_arr = crds.req_postions(1:(end-1),1)';
                        crds.y_arr = crds.req_postions(1:(end-1),2)';
                        exp_data.coords = [crds.x_arr' , crds.y_arr'];
                    else
                        error('Unregonized mat file data structure.');
                    end
                else
                    if strcmpi(crds_EXT,'.hdf5')
                        % as from U:\Ptycho\PtychoShelves_EM-1.0\ptycho\+scans\+positions\hdf5_pos.m
                        ps = h5read(crds_file_name,'/probe_positions_0');
                        ppX = ps(:,1)*1e-10; % Angstrom to meter --> ppX is read in column 1
                        ppY = ps(:,2)*1e-10; % --> ppY is read in column 2

                        crds.x_arr  = -ppY(:).'; % Note column 1 is -1*read in column 2
                        crds.y_arr  = -ppX(:).'; % Note column 2 is -1*read in column 1
                        
                    else
                        error('Unregonized file format for scan coordinates');
                    end
                end
            end
            clear dp_set
        else
            exp_data = loadHPL(loadup);
            if ~isfield(exp_data,'diffPats')
                if isfield(exp_data,'diffraction')
                    if isfield(exp_data.diffraction,'micrograph')
                        exp_data.diffPats = exp_data.diffraction.micrograph; % when loaded has something like:     micrograph: [2025×512×512 single]
                        exp_data.diffraction = rmfield(exp_data.diffraction,'micrograph');
                        permute_order_of_loaded_dps = [2, 3, 1];
                        exp_data.diffPats = permute(exp_data.diffPats, permute_order_of_loaded_dps);    % this makes it look like as seen in azorus, but with index as last element.
                    else
                        error('Cannot find diffraction data in loaded HP file');
                    end                       
                end
            end
        end
        if any(strcmpi(set_params.Properties.VariableNames, 'ignore_pixels')) 
            if iscell(set_params.ignore_pixels) && cell_is_charorstring(set_params.ignore_pixels{set_number})
                [basef, fn, fext] = fileparts(loadup);
                if strcmpi(set_params.ignore_pixels{set_number},'_MASK')
                    maskname = fullfile(basef,[fn,'_MASK',fext]);
                else
                    maskname = fullfile(basef,set_params.ignore_pixels{set_number});
                end
                mask_in = load(maskname);
                if isfield(mask_in,'mask')
                    exp_data.mask = mask_in.mask;
                else
                    loadedfnames = fieldnames(mask_in);
                    exp_data.mask = mask_in.(loadedfnames{1});
                    warning('Loaded first field of mask file, this might not be the correct one.');
                end
                clear mask_in;
                if ~isempty(permute_order_of_loaded_dps)
                    exp_data.mask = permute(exp_data.mask, permute_order_of_loaded_dps);  % make sure this the same as above.
                end
            end
        end
    else
        exp_data = loadHPL(fullfile(experi_def.localsummary_dir ,'mini_test_set.HPL'));
        exp_data.meta.angular_ps = 0.5E-6;
    end
    
    if exist('trans_coords','var')
        if ~exist('crds','var')
            crds.x_arr = exp_data.coords(:,1)';
            crds.y_arr = exp_data.coords(:,2)';
        end
        switch(trans_coords)
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
                error('Invalid trans_coords');
        end
    end
    
    
    if  any(strcmpi(set_params.Properties.VariableNames, 'simulated_dose')) && ...
        contains_num(set_params.simulated_dose(set_number))
        fprintf('Adjusting Dose to %d e/A^2...\n',set_params.simulated_dose(set_number));
        % Synthetically adjust the diffraction patterns to have the dose:
        % First get native average dose:
        % - Average num of electrons per dp:
        ave_counts_per_dp = sum(exp_data.diffPats(:))/size(exp_data.diffPats,3);
        % assume regular square grid on average:
        % (In recent simulations eg,
        % D:\Local\GlobusLandingArea\Python20200530\AuCubes07_TDS01.py,
        % have been doing something along the lines of:
        %{
        final_dpSize = 256
        recon_pixel_size = 0.018 #nm
        Diameter_Aim = 10*(recon_pixel_size*final_dpSize)/3; # Angstrom
        diam_frac = 0.17;
        step_size = diam_frac * Diameter_Aim
        %}
        patch_area = set_params.step_size_A(set_number)^2;
        initial_dose = ave_counts_per_dp/patch_area;
        exp_data.diffPats = (set_params.simulated_dose(set_number)/initial_dose) * exp_data.diffPats;
        exp_data.diffPats = add_ele_noise(exp_data.diffPats); % adds poissionian noise.
        % The above is of course only correct if the initial set was noise
        % free...
        fprintf('Dose Adjusted.\n');
%         figure; imagesc(exp_data.diffPats(:,:,1)); axis image
    end
        
    
    if any(strcmpi(set_params.Properties.VariableNames, 'adjacent_params')) && ...
             contains_num(set_params.adjacent_params(set_number))
        set_number = set_number + set_params.adjacent_params(set_number); %#ok<FXSET>
        % Now that's what I call potentially confusing and ugly
        % programming... changing the loop counter inside the loop. 
    end
    if ~experi_def.summary_table.use_spreadsheet_experi_params
    % use built in calib
        scale_posn =1;
        coords_nm = exp_data.coords_nm';
        collection_params = exp_data.meta;
        stored_recipe = exp_data.pi_recipe.recipe;        
        DP_SIZE = exp_data.meta.detector_shape(1);
        BeamVoltage = collection_params.highVoltage*1e3;
        disp('Im here');
        rotation_angle = 180*(collection_params.affineRotation/pi);
    else
    % or use from spreadsheet:    
        % coords_nm = transpose(exp_data.coords * [[0 -1];[1 0]] * nm_per_V);
        % figure; plot(coords_nm(1,:), coords_nm(2,:));
        % hold on; plot(coords_nm(1,1:10), coords_nm(2,1:10),'rx');
        collection_params = [];
        stored_recipe = check_pty_recipe(0);
        BeamVoltage = set_params.kV(set_number)*1e3;
        rotation_angle = set_params.meas_rot_est(set_number);
        beam_half_angle = set_params.beam_half_angle_mrad(set_number)*1e-3;
        illum_diam = set_params.illum_diam(set_number); %nm
        recon_pix_size_nm = set_params.recon_pix_size_pm(set_number)/1000; 
        nm_per_V = set_params.nm_V(set_number);
        scale_posn = set_params.scale(set_number);
        coords_nm = transpose(exp_data.coords * [[0 -1];[1 0]] * nm_per_V); 
        % swaps row 1 and 2, and multiples resulting first by -1. Azorus
        % for example works top down, with first row of coords being the
        % row number, as used in conventional image space.
%         if isfield(set_params,'set_number') && isnumeric(set_params.DP_SIZE(set_number))
%             DP_SIZE = set_params.DP_SIZE(set_number);
%         else
%             DP_SIZE = 256;
%         end
        DP_SIZE = size(exp_data.diffPats,1);
    end
    beam = RelatBeam(BeamVoltage, 0);
    lambda = beam.BeamLambdaRel;
    % prep the main data:
    % should make load HPL not do initial permute? -> Done March 23:
    % exp_data.diffPats = permute(exp_data.diffPats,[3 2 1]);
    % use built in params:
    if ~experi_def.summary_table.use_spreadsheet_experi_params
    % use built in calib 
        angle_per_pix = exp_data.meta.angular_ps;
        recon_pix_size = lambda / ( DP_SIZE*angle_per_pix);
        recon_pix_size_nm = recon_pix_size/1e-9;
    else
        % use known pixel size
        % recon_pix_size_nm = 18e-3;
        recon_pix_size = recon_pix_size_nm*1e-9; % nm
        angle_per_pix = lambda / ( DP_SIZE*recon_pix_size);
    end
    
    
    try
%         parfor (recipe_var = experi_def.recipe_table.use_variations)    
        for recipe_var = experi_def.recipe_table.use_variations    
    %       with parfor we need initialize temporary variables within the loop
    %       to be totally clear as to what's going on and get no differences
    %       between for and parfor loops:

            drift_angle_deg =0;
            drift_mag_nm_during_acquire =0;
            det_subregion =[];
            prime_recon_fname = experi_def.recipe_table.prime_probe_num; %TO DO: fix this to be give a sensible name....
            if exec_params.parallel
                last_saved_file = []; % In parfor which file is the last saved is unclear: can't use 'step_previous' in parfor calls. TO DO: graceful way of dealing with this in code below.
                last_montage = []; % There is no last montage in parallel sitation... For that matter plotting anything in parfor rapidly becomes an issue....            
            end
            lsqml_ver = 1; % by defaults goes to old version.
%             recipe_calc();
%                 feval(fcn);
    
    
%     function recipe_calc()
    %   %% Setup 
        PS_recon = ptychoset_mult;
        save_dir = base_dir;
        data_setname = 'localcache';
        if exec_params.use_pi_box
            PS_recon.setnames(save_dir, data_setname);
            PS_recon.MakeDirs;  
        end

        if cell_is_true(set_params.use_drift(set_number))
           if set_params.calc_drift(set_number)
               before_image = fullfile(base_dir, set_params.before_drift_im{set_number});
               after_image  = fullfile(base_dir, set_params.after_drift_im{set_number});
%                [drift_mag_nm_between_bef_and_after, drift_angle_deg, drift_mag_nmpersec] = ptychoset_mult.calc_drift(before_image, after_image);
               [~, drift_angle_deg, drift_mag_nmpersec] = ptychoset_mult.calc_drift(before_image, after_image);
               drift_mag_nm_during_acquire = drift_mag_nmpersec*(set_params.exposure_time(set_number)/1000)*length(coords_nm);
           else
               drift_mag_nm_during_acquire = set_params.drift_mag_nm(set_number);
               drift_angle_deg = set_params.drift_angle_deg(set_number);
           end
        end

        ini_recon_params = PS_recon.setup_recon_params('struct',struct());
        ini_recon_params.dpSize = DP_SIZE;
        ini_recon_params.angle_per_pix = angle_per_pix;
        ini_recon_params.BeamVoltage = BeamVoltage;
        ini_recon_params.calc_probe = 0;
        PS_recon.setup_recon_params('struct',ini_recon_params);
        PS_recon.recon_dat.p.crds_file_name = crds_file_name;

        % Set the Positions up for passing to pi-box
        PS_recon.dps.scan_posns = [coords_nm(1,:)*scale_posn;...
                                   coords_nm(2,:)*scale_posn];
        PS_recon.dps.crds_file_name = crds_file_name;
        if isfield(exp_data,'attrs') && isfield(exp_data.attrs,'meshParams') && isfield(exp_data.attrs.meshParams,'shape')
            PS_recon.dps.attrs.meshParams.shape = exp_data.attrs.meshParams.shape;
        elseif iscell(set_params.scan_shape(set_number)) && ischar(set_params.scan_shape{set_number})
            PS_recon.dps.attrs.meshParams.shape = str2num(set_params.scan_shape{set_number});
        else
            % could port the code I wrote in python to find out the scan shape
            % from the positions, assuming a raster pattern layout...
        end

        % Was any rotation applied to the scan coords in the first
        % place? Note the that drift measurement would have been taken
        % in a non rotated frame.                       
%         if exist('collection_params','var') && isfield(collection_params, 'affineRotation')
        if ~(experi_def.summary_table.use_spreadsheet_experi_params) && isfield(collection_params, 'affineRotation')
            preapplied_rotation_angle = 180*(collection_params.affineRotation/pi);
        elseif contains_num(set_params.rot_applied(set_number))
            preapplied_rotation_angle = set_params.rot_applied(set_number);
        else
            preapplied_rotation_angle = 0;
        end

        ini_recon_params.rot_angle_deg = -1*preapplied_rotation_angle;
        PS_recon.setup_recon_params('struct',ini_recon_params);
        PS_recon.rotate_dp_posns('scan_posns'); % this rotates scan positiions and stores them in '.rotated_posns'
         if islogical(subregion) && subregion
             if isempty(set_params.central_RCsize(set_number))
                 [~, det_subregion] = PS_recon.create_subset;
%                  [sub_region_pnts, det_subregion] = PS_recon.create_subset;
             else
                 if iscell(set_params.central_RCsize) && ischar(set_params.central_RCsize{set_number})
                     RCsize = str2num(set_params.central_RCsize{set_number});
%                      [sub_region_pnts, det_subregion] = PS_recon.create_subset(RCsize);
                     [~, det_subregion] = PS_recon.create_subset(RCsize);
                 else
                     error(['Cannot determine the subregion size from input']); %,\n',...
    %                           'Reverting to manual selection through GUI.\n']);
    %                  [sub_region_pnts, subregion] = PS_recon.create_subset;
                 end

             end
         end

         % might want to break around here to save points and DPs in a file, so
         % that we don't have to load the superset each time in order to get
         % the subset, e.g. for first implementation of PtychoShelvesSet:
         %{


         %}

        if cell_is_true(set_params.use_drift(set_number))
    %     if ~isempty(set_params.use_drift(set_number)) && logical(set_params.use_drift(set_number)) 
             % Apply drift corrections, rotating back to non-rotated frame.
    %              drift_mag_nm_during_acquire is the drift of the LAST diffraction exposure
    %              position.

             PS_recon.drift_scan_posns(drift_mag_nm_during_acquire, drift_angle_deg, 'rotated_posns'); % this stores results in 'drifted_posns'
             % now we want rotate back to the pre-applied position:
             ini_recon_params.rot_angle_deg = preapplied_rotation_angle;
             PS_recon.setup_recon_params('struct',ini_recon_params);
             PS_recon.rotate_dp_posns('drifted_posns'); % this rotates scan positiions and stores them in '.rotated_posns'
        else

             ini_recon_params.rot_angle_deg = preapplied_rotation_angle;
             PS_recon.setup_recon_params('struct',ini_recon_params);
             PS_recon.rotate_dp_posns('rotated_posns');
             % there's a chance we just went back and forward by 0 degrees, but
             % at least the results are consistently in rotated posns now.
        end

        if ~isempty(rotation_angle)
            % apply a rotation to the positions
            ini_recon_params.rot_angle_deg = rotation_angle - preapplied_rotation_angle;
        else
            ini_recon_params.rot_angle_deg = 0;
        end
        PS_recon.setup_recon_params('struct',ini_recon_params);
        PS_recon.rotate_dp_posns('rotated_posns'); % this rotates scan positiions and stores them in '.rotated_posns'  

        if ~isempty(subregion) && ((subregion ~= 0) || (det_subregion ~=0))
            PS_recon.dps.rotated_posns = PS_recon.dps.rotated_posns(:,det_subregion);
            PS_recon.dps.scan_posns    = PS_recon.dps.scan_posns(:,det_subregion); 
        end

        PS_recon.create_RC_coords('rotated_posns'); % this converts '.rotated_posns' to RowColumn coords for input into pi-box and creates the default object_0
        default_object = PS_recon.recon_dat.object_0;
   
%%
        PS_recon.recon_dat = [];
        PS_recon.setup_recon_params('struct',ini_recon_params);
        
        % this needs to be here, as we're now manipulatinng the diff pats
        % before passsing to the lsqml recipes.
        if ~isempty(subregion) && (subregion ~= 0)
            PS_recon.dps.ims_proc = exp_data.diffPats(:,:,det_subregion); % this would go wrong if expdata has data in wrong order....
            PS_recon.dps.masks = exp_data.mask(:,:, det_subregion);
        else
            PS_recon.dps.ims_proc = exp_data.diffPats;
            if isfield(PS_recon.dps,'masks')
                PS_recon.dps.masks = exp_data.mask;
            end
        end 
        
%         %%        
        new_recipe = stored_recipe;
        stored_param_names = fieldnames(stored_recipe.params(1));
        recipe_variations_param_names = recipe_variations.Properties.VariableNames;

        if ~isempty(recipe_variations.sub_steps{recipe_var})
            n_sub_steps = length(recipe_variations.sub_steps{recipe_var});
        else
            n_sub_steps = 1;
        end
        
        for sub_step_i = 1:n_sub_steps % for setting up when using recipe with substeps listed e.g., '4, 5, 6'

            for recipe_param_i = 1:length(recipe_variations_param_names)
                variation_field_name = recipe_variations_param_names{recipe_param_i};
                if iscell(recipe_variations.(variation_field_name))
                    variation_value = recipe_variations.(variation_field_name){recipe_var + sub_step_i -1};
                else
                    variation_value = recipe_variations.(variation_field_name)(recipe_var + sub_step_i -1);
                end
                
                if ~isempty(variation_value) && ~isnumeric(variation_value) && ~islogical(variation_value)
                    var_value = str2num(variation_value); %#ok<ST2NM>
                    if isempty(var_value)
                        var_value = variation_value;
                    end
                    variation_value = var_value;
                end
                
                if ~isempty(variation_value)
%                     if isfield(stored_recipe.params(1), variation_field_name)
                        new_recipe.params(sub_step_i).(variation_field_name) = variation_value;                    
%                     end
                    % if the field is empty leaves it as it was in loaded
                    % recipe
                end
            end
            % and the number of iterations to use: (oriignallly below were
            % cell refs {} on iterations...
            if ~isempty(recipe_variations.iterations(recipe_var + sub_step_i -1))
                if ~isnumeric(recipe_variations.iterations(recipe_var + sub_step_i -1))
                    iterations = str2num(recipe_variations.iterations(recipe_var + sub_step_i -1)); %#ok<ST2NM>
                else
                    iterations = recipe_variations.iterations(recipe_var + sub_step_i -1);
                end
                new_recipe.iterations(sub_step_i) = iterations;
            end
        
        end
        
        % Delete the later steps in the recipe that are not to be used to
        % form recipes:
        if n_sub_steps < length(new_recipe.params)
            new_recipe.params((n_sub_steps+1):end) = [];
            new_recipe.iterations((n_sub_steps+1):end) = [];
        end
        
                
        % Now sort out what to do with the probe and object guesses:
        % place fields within new_recipe to say whether to use probe or object
        % from a previous iteration or model:
        new_recipe.probe_model = recipe_variations.probe_model{recipe_var};
        prb_loaded_from_earlier_step = false;
        switch new_recipe.probe_model
            case 'step_prime'
                primary_dat = load([prime_recon_fname,'.mat']);
                PS_recon.recon_dat.probe_pf = primary_dat.recon_dat.probe_int{end}; % currently need to transpose this when using pi-box: TODO sort this shit out.
                fprintf('Using probe model from %s for computation.\n',prime_recon_fname);
                new_recipe.probe_model_file = prime_recon_fname;
                prb_loaded_from_earlier_step = true;
            case 'step_previous'
                last_dat = load([last_saved_file,'.mat']);
                PS_recon.recon_dat.probe_pf = last_dat.recon_dat.probe_int{end}; % currently need to transpose this when using pi-box: TODO sort this shit out.
                fprintf('Using probe model from %s for computation in step.\n',last_saved_file);
                new_recipe.probe_model_file = last_saved_file;
                prb_loaded_from_earlier_step = true;
            case {'step_numbered','step_numbered_masked'}
%                 reconsavedpath = base_dir;
                if ischar(set_params.probe_load_name{recipe_var})
                    base_load_name = set_params.probe_load_name{recipe_var}; 
%                     base_load_name = set_params.save_name{recipe_variations.probe_expt_num(recipe_var)};
                else
                    base_load_name = save_name;
                end
                num_savename = sprintf('%s_recipe%03.0f',base_load_name, recipe_variations.probe_recipe_num(recipe_var));
                
                % TO DO need a programatic way to implement the below, to allow
                % probe to be pulled in from a different experiment number, which might thus have a different store location to the current savepath.
                % The below would be / is a hazard as it partly defeats the
                % idea of the forces_same_expt_number switch,
                
%                 probe_location_path = set_params.recon_savepath{recipe_variations.probe_expt_num(recipe_var)};
%                 destination_file = fullfile(probe_location_path, num_savename);
                % NOTE: IF PROBE IS LOADED FROM A DIFFERENT EXPERIMENT SET, NO
                % CHECKS OR SCALING IS DONE HERE TO GET THE PIXEL SIZE
                % RIGHT, AND SCALE APPROPRIATELY. THUS USE WITH CAUTION!
                
                destination_file = fullfile(reconsavepath, num_savename);
                numbered_name = [destination_file,'.mat'];
                num_dat = load(numbered_name);
                
                
                PS_recon.recon_dat.probe_pf = num_dat.recon_dat.probe_int{end}; % currently need to transpose this when using pi-box: TODO sort this shit out.
                fprintf('Using probe model from %s for computation in step.\n',numbered_name);
                new_recipe.probe_model_file = numbered_name;
                prb_loaded_from_earlier_step = true;
                % FOR LSQML need to make sure that probe data is in correct
                % field name... could also not do the below, and use
                % function later to write probe (loaded above) to a known
                % location, where the data is in that file with the correct
                % fieldname.
%                     PS_recon.recon_dat.ptychoshelves.initial_probe_file = numbered_name;
                if isfield(new_recipe.params,'Probe_Init_Recentre') && ...
                    ~isnan(new_recipe.params.Probe_Init_Recentre) && ...
                    new_recipe.params.Probe_Init_Recentre > 0 

                    ixs = [{':',':'},num2cell(ones(1,ndims(PS_recon.recon_dat.probe_pf)-2))];
                    prb = PS_recon.recon_dat.probe_pf(ixs{:});

                    centre_of_prb = FindCentre(abs(prb));
                    PS_recon.recon_dat.probe_pf = Recentre(prb, centre_of_prb);                
                end


                if isfield (new_recipe.params,'Probe_Init_Circ_Mask') && ...
                    ~isnan(new_recipe.params.Probe_Init_Circ_Mask) && ...
                     new_recipe.params.Probe_Init_Circ_Mask > 0
%                 if strcmpi(new_recipe.probe_model,'step_numbered_masked')
                    prbsize = size(PS_recon.recon_dat.probe_pf);
                    [px, py] = meshgrid( (1:prbsize(2)) - round(prbsize(2)/2) - 1 , ...
                                         (1:prbsize(1)) - round(prbsize(1)/2) - 1);
                    probe_mask = circ(px, py, new_recipe.params.Probe_Init_Circ_Mask*min(prbsize));
                    PS_recon.recon_dat.probe_pf = PS_recon.recon_dat.probe_pf.*probe_mask;
                end

            case 'builtin'
                fprintf('Using probe guess already computed and stored in hpl or mat file.\n');
                if isfield(exp_data,'probeModel')
                    PS_recon.recon_dat.probe_pf = exp_data.probeModel;
                elseif isfield(exp_data,'probe0')
                    PS_recon.recon_dat.probe_pf = exp_data.probe0;
                else
                    error('Cannot find field in data containing probe model');
                end
            case 'xlcompute'
                fprintf('Computing probe guess from values in excel file.\n');
                if isempty(set_params.probe_mod_defocus(set_number))
                    defocus = (set_params.under_or_over(set_number))*((illum_diam/2)/beam_half_angle)*1e-9;
                else
                    defocus = set_params.under_or_over(set_number) * set_params.probe_mod_defocus(set_number) * 1e-9;
                end
                % in this file set_params is a table, so fieldnames will
                % give incorrect answer. I kind of wish isfield call on
                % a table threw an error or a warning... maybe make
                % ammended isfield function or local override?
                if any(strcmpi(fieldnames(set_params), 'Cs'))  && ~isempty(set_params.Cs(set_number)) 
                    Cs = 1e-3*set_params.Cs(set_number); % 1e-3 as it is given in mm in spreadsheet.
                else
                    Cs = [];
                end
                
                % For record keeping put parameters here:
                PS_recon.recon_dat.beam_half_angle = beam_half_angle;
                PS_recon.recon_dat.defocus = defocus;
                PS_recon.recon_dat.Cs = Cs;                                
                
                % Use perfect probe with ptychoset:                                
                PS_recon.recon_dat.probe_pf = ...
                     ptychoset_mult.perfect_probe(DP_SIZE, recon_pix_size, defocus, beam_half_angle, lambda, [], Cs);
            case 'xlcompute_powernorm'
                fprintf('Computing probe guess from values in excel file and doing power norm.\n');
                if isempty(set_params.probe_mod_defocus(set_number))
                    defocus = (set_params.under_or_over(set_number))*((illum_diam/2)/beam_half_angle)*1e-9;
                else
                    defocus = set_params.under_or_over(set_number) * set_params.probe_mod_defocus(set_number) * 1e-9;
                end
                if any(strcmpi(fieldnames(set_params), 'Cs'))  && ~isempty(set_params.Cs(set_number)) && ~isnan(set_params.Cs(set_number)) 
                    Cs = 1e-3*set_params.Cs(set_number); %1e-3 as it is given in mm in spreadsheet.
                else
                    Cs = [];
                end
                % For record keeping put parameters here:
                PS_recon.recon_dat.beam_half_angle = beam_half_angle;
                PS_recon.recon_dat.defocus = defocus;
                PS_recon.recon_dat.Cs = Cs;   
                
                
                % Use perfect probe with ptychoset:
                PS_recon.recon_dat.probe_pf = ...
                     ptychoset_mult.perfect_probe(DP_SIZE, recon_pix_size, defocus, beam_half_angle, lambda, [], Cs);
                
                average_root_dp_intens = sum(sqrt(single(PS_recon.dps.ims_proc(:))))/size(PS_recon.dps.ims_proc,3);
                % do our own scaling out of interest to see if this helps..
                model_probe_intens = sum(abs(PS_recon.recon_dat.probe_pf(:)));
                model_probe_dp_intens = model_probe_intens*DP_SIZE^2;
                model_intens_factor = average_root_dp_intens /model_probe_dp_intens;
                PS_recon.recon_dat.probe_pf = PS_recon.recon_dat.probe_pf * model_intens_factor;
                
            otherwise
                fprintf(['No or unknown probe specification given.\n',...
                         'Using default built into model.\n']);
                PS_recon.recon_dat.probe_pf = exp_data.probeModel;    
        end
        
        if isfield(new_recipe,'scale_probe_pix') && ...
                ~isempty(new_recipe.scale_probe_pix) && ...
                 new_recipe.scale_probe_pix > 0 && ...
                 prb_loaded_from_earlier_step
             
             % e.g. if pixel size of probe model supplied is 13 pm, and
             % that required for the current model is 17 pm, then in this
             % case FOV of supplied probe is 13 nm (1000 px), but our
             % current probe has FOV of 17 nm. Thus we need would need to
             % scale probe image by (13/17) then padd with some value on
             % outside. Alternatively if other way around wuuld need to 
             % crop back to the same size.
             PS_recon.recon_dat.probe_pf = ...
                    imresize(real(PS_recon.recon_dat.probe_pf), new_recipe.scale_probe_pix) + ...
                 1i*imresize(imag(PS_recon.recon_dat.probe_pf), new_recipe.scale_probe_pix);
             PS_recon.recon_dat.probe_pf = utils.crop_pad( PS_recon.recon_dat.probe_pf, ...
                 size(PS_recon.recon_dat.probe_pf), 0);
             
        end

%         forgot to comment this out before!
%         if exec_params.use_pi_box && prb_loaded_from_earlier_step
%             PS_recon.recon_dat.probe_pf = PS_recon.recon_dat.probe_pf.';  %TODO: Get to the bottom of this hack.
%         end
%%        

        new_recipe.object_model = recipe_variations.object_model{recipe_var};
        obj_loaded_from_earlier_step = false;
        
        % defaults for ptychoset:
        PS_recon.recon_dat.ptychoshelves.model_object = false; % causes a load from file below.
        PS_recon.recon_dat.ptychoshelves.object_type = 'rand';        
        PS_recon.recon_dat.ptychoshelves.initial_iterate_object_file{1} = '';
        PS_recon.recon_dat.ptychoshelves.distrib_objphase = false;
        PS_recon.recon_dat.ptychoshelves.mean_filter_obj = [];
        PS_recon.recon_dat.ptychoshelves.tie_unwrap_obj = false;
%         isfield(recipe_variations,'compute_cumu_phases') % don't do this
%         as recipe_variations is a table
        if any(strcmpi(fieldnames(recipe_variations), 'compute_cumu_phases')) && ...
            ~isnan(recipe_variations.compute_cumu_phases(recipe_var))
            PS_recon.recon_dat.ptychoshelves.compute_cumu_phases = ...
               recipe_variations.compute_cumu_phases(recipe_var);
        else
            PS_recon.recon_dat.ptychoshelves.compute_cumu_phases = false;
        end
        PS_recon.recon_dat.ptychoshelves.cumu_phase_regulation = ...
            recipe_variations.cumu_phase_regulation(recipe_var);
        
        switch new_recipe.object_model
            case 'step_prime'
                primary_dat = load([prime_recon_fname,'.mat']);
%                 PS_recon.recon_dat.object_0 = primary_dat.recon_dat.obj_int{end}.'; % currently need this when using pi-box: TODO sort this shit out.
                PS_recon.recon_dat.object_0 = primary_dat.recon_dat.obj_int{end};
                new_recipe.object_model_file = prime_recon_fname;
                obj_loaded_from_earlier_step = true;
            case 'step_previous'
                last_dat = load([last_saved_file,'.mat']);
%                 PS_recon.recon_dat.object_0 =
%                 last_dat.recon_dat.obj_int{end}.'; % currently need this when using pi-box: TODO sort this shit out.
                PS_recon.recon_dat.object_0 = last_dat.recon_dat.obj_int{end};
                new_recipe.object_model_file = last_saved_file;
                obj_loaded_from_earlier_step = true;
            case 'step_numbered'
%                 reconsavedpath = base_dir;
                if ischar(set_params.object_load_name{recipe_var})
                    base_load_name = set_params.object_load_name{recipe_var}; %% good!
                else
                    base_load_name = save_name;
                end
                num_savename = sprintf('%s_recipe%03.0f',base_load_name, recipe_variations.object_recipe_num(recipe_var));
                destination_file = fullfile(reconsavepath, num_savename);
                numbered_name = [destination_file,'.mat'];
                num_dat = load(numbered_name);
%                 PS_recon.recon_dat.object_0 = num_dat.recon_dat.obj_int{end}.';% currently need this when using pi-box: TODO sort this shit out.
                PS_recon.recon_dat.object_0 = num_dat.recon_dat.obj_int{end};
                fprintf('Using Object model from %s for computation in step.\n',numbered_name);
                new_recipe.object_model_file = numbered_name;
                obj_loaded_from_earlier_step = true;
                % for LSQML 2
                PS_recon.recon_dat.ptychoshelves.model_object = false;
                PS_recon.recon_dat.ptychoshelves.object_type = [];
                % TO DO ** make sure object is saved in correct fieldname within the mat file.                 
                PS_recon.recon_dat.ptychoshelves.initial_iterate_object_file{1} = numbered_name;
                PS_recon.recon_dat.ptychoshelves.distrib_objphase = recipe_variations.distrib_objphase(recipe_var);                
                PS_recon.recon_dat.ptychoshelves.mean_filter_obj = recipe_variations.mean_filter_obj(recipe_var);
%                 PS_recon.recon_dat.ptychoshelves.mean_filter_obj = str2num(PS_recon.recon_dat.ptychoshelves.mean_filter_obj{:}); %#ok<ST2NM>
                PS_recon.recon_dat.ptychoshelves.tie_unwrap_obj  = recipe_variations.tie_unwrap_obj(recipe_var);
            case 'blank'
                % use default.
                PS_recon.recon_dat.object_0 = default_object;
            case 'flat'
                % this is one the models from the ptychoshelves:
                PS_recon.recon_dat.ptychoshelves.model_object = true;
                PS_recon.recon_dat.ptychoshelves.object_type = 'flat';
                PS_recon.recon_dat.ptychoshelves.initial_iterate_object_file{1} = '';
                PS_recon.recon_dat.object_0 = [];
            case 'flat0'
                % this is one the models from the ptychoshelves:
                PS_recon.recon_dat.ptychoshelves.model_object = true;
                PS_recon.recon_dat.ptychoshelves.object_type = 'flat0';
                PS_recon.recon_dat.ptychoshelves.initial_iterate_object_file{1} = '';
                PS_recon.recon_dat.object_0 = []; 
            case 'rand'
                % this is one the models from the ptychoshelves:
                PS_recon.recon_dat.ptychoshelves.model_object = true;
                PS_recon.recon_dat.ptychoshelves.object_type = 'rand';
                PS_recon.recon_dat.ptychoshelves.initial_iterate_object_file{1} = '';
                PS_recon.recon_dat.object_0 = [];
                
            otherwise
                fprintf(['Unknown object specification given.\n',...
                         'Using default blank object starting model.\n']);
                PS_recon.recon_dat.object_0 = default_object;
        end
        
%         if exec_params.use_pi_box && obj_loaded_from_earlier_step
%             PS_recon.recon_dat.object_0 = PS_recon.recon_dat.object_0.';
%         end        
        
        new_recipe.circular_mask_rad = recipe_variations.circular_mask_rad(recipe_var); 
        if new_recipe.circular_mask_rad > 0
            DPmask = circularMaskGen(DP_SIZE, new_recipe.circular_mask_rad);
        else
            DPmask = [];
        end
                 
%%
        % Sanity check and prepare recipe:
        fprintf('Executing Expt %.0f, Recipe %.0f ...\n', set_number, recipe_var);
        fprintf('Recon Pixel Size %3.1f pm \n',PS_recon.recon_dat.recon_pixel_size/1e-12);
        n_steps = length(new_recipe.params);
        t1 = now;
        
        
        %% THE RECIPE EXECUTION:
        
        show_newposns = 0;
        
        try
        
            if exec_params.use_lsmql                
                iter = 1; % TO DO:
                % Need to figure out how best sequences of params can be used in lsmql code            
                % *** Note that the order of the operations below matters*** , and I
                % don't know exactly what the pi-box is doing at the moment...
                % Note: Should put this above somewhere so it can deal with
                % subregions properly
                % Update ... June 2021 - Ptychoshelves, i.e. lsqml_ver = 2,
                % should make seequences much easier now!
                if isfield(new_recipe.params,'code') && cell_is_charorstring(new_recipe.params(iter).code)
                    switch lower(new_recipe.params(iter).code)
                        case 'ps_em'
                            lsqml_ver = 2;
                        otherwise
                            lsqml_ver = 1;
                    end                    
                end

                if lsqml_ver == 1
                    fprintf('lsqml_ver == 1 used to implement an earlier version of LSQML in ptychoshelves. See archives.\n');
                    
                elseif lsqml_ver == 2
                    local_params = exec_params;
                    local_params.recon_savename = sprintf('%s_recipe%03.0f',save_name, recipe_var);
                    probe_savename = fullfile(set_params.recon_savepath{set_number},sprintf('probe_expt%3.0f.mat',set_number));
                    PS_recon.save_probe_pf(probe_savename); % this saves the probe generated by ptychoset into a location, and puts that location path into PS_recon.recon_dat.ptychoshelves.initial_probe_file 
                    PS_recon.recon_dat.ptychoshelves.crds_file = crds_file_name;
                    [p, eng] = PS_recon.conv_ptychoset_cSAXS_MSlc(set_params, set_number, new_recipe, iter, local_params); % set_params.trans_coords(set_number)                                
                    [p, ~] = core.append_engine(p, eng); 
                    outputs = core.ptycho_recons(p);

                    PS_recon.recon_dat.probe_int{1} = outputs.probes;
                    PS_recon.recon_dat.obj_int{1} = outputs.object{1};
                    PS_recon.recon_dat.recon_param = p;
                    PS_recon.recon_dat.error_metric = outputs.error_metric;
                    PS_recon.recon_dat.reconstruction_engine_time = outputs.computation_time;

                    % Create mimimal version of output from solver:
                    engines = rmfield(outputs.engines{1},{'object_final','probes_final'});
                    outp = rmfield(outputs,{'fmag','probes','object'});
                    outp.engines = engines;
                    PS_recon.recon_dat.shelves_op = outp; % add this in so that it is saved in mat file
                end           



                % At one point I wondered if oddities were happening in PSI
                % functions due to use of persistent global variables in lsqml
                % code, so I blew the cobwebs away here. I did not see any issues, 
                % but in case some oddity happens again uncomment this... Watch
                % for deadly global var name clashes etc...
    %             clear functions %#ok<CLFUNC>
    %             clear core.ptycho_solver
    %             clear core
    %             clear global
    %             clear mex %#ok<CLMEX>
            end


            t2 = now; % timer for recipe execution

            PS_recon.recon_dat.exec_start = t1;
            PS_recon.recon_dat.exec_time = (t2 - t1)*(24*60*60);
            % should make a function to do nice contrast setting on data:
            % caxis(absObjAxe,[1 1.5]);
            % caxis(angObjAxe,[-0.5 0.2]);

            recon_dat = PS_recon.recon_dat;
            % add in the recipe data:
            recon_dat.recipe = new_recipe;
            recon_dat.rotation_angle = rotation_angle;
            % add a field in with the propogated image:
            obj_size = size(recon_dat.obj_int{end});
            if length(obj_size) > 2
                nslices = obj_size(end);
            else
                nslices = 1;
            end
            if nslices > 1
                obj_inds = num2cell(ones(1,length(obj_size)));
                obj_inds(1:2) = {':'};
                im0 = squeeze(recon_dat.obj_int{end}(obj_inds{:}));
                for slc = 2:obj_size(end)
                    im1_illum = slice_progate(im0, ...
                        recon_dat.recon_pixel_size,...
                        recon_dat.lambda, ...
                        recon_dat.recipe.params.delta_z_angstrom*1e-10);
                    obj_inds{end} = slc;
                    im0 = im1_illum .* squeeze(recon_dat.obj_int{end}(obj_inds{:}));
                end
                recon_dat.propogated_obj = im0;
            end

            recon_savename = sprintf('%s_recipe%03.0f',save_name, recipe_var);
            destination_file = fullfile(reconsavepath, recon_savename);
            if ~exist(destination_file, 'file')
                save([destination_file,'.mat'],'recon_dat','-v7.3');
            else            
                new_dest_file = [destination_file, '_',datestr(now,30)]; 
                fprintf(['Reconstruction savename already exists.\n',...
                    'Stored in new file name: %s.mat'],destination_file);
                save([new_dest_file,'.mat'],'recon_dat','-v7.3');
            end
            last_saved_file = destination_file;                 
            tif_summary_savename = [last_saved_file,'_summary.tif'];
            
            text_description = recon_savename;
            % text_description = SerializeReconInfo(PS_recon.recon_dat));
           
            
            if nslices == 1            
                savesummaryimage3(tif_summary_savename , ...
                         recon_dat.obj_int{end},...
                         recon_dat.probe_int{end},...
                         -recon_dat.rotation_angle,...
                         text_description);
            else
%                 savesummaryimage3(tif_summary_savename , ...
%                          recon_dat.propogated_obj,...
%                          recon_dat.probe_int{end},...
%                          -recon_dat.rotation_angle,...
%                          text_description);                
                savesummaryimage3(tif_summary_savename, ...
                  recon_dat.obj_int{end}, recon_dat.probe_int{end}, ...
                 -recon_dat.rotation_angle,...
                 text_description, struct('phase_pres_method','propa',...
                       'recon_pix_size', recon_dat.recon_pixel_size, ...
                       'lambda', recon_dat.lambda ,...
                       'Dz', recon_dat.recipe.params.delta_z_angstrom*1e-10));
            end
    
            if exec_params.copy_to_monitor_dir     
                try     
                    copyfile(tif_summary_savename, exec_params.local_monitor_dir);
                catch
                    fprintf('Something went wrong with copy ... No biggie\n');
                end
            end
            recent_im = imread(tif_summary_savename);
            if exec_params.show_last_recon && ~exec_params.parallel
                if  ~exist('last_montage','var') || ~(ishandle(last_montage) && findobj(last_montage,'type','axes')==last_montage)
                    f1 = figure;
                    last_montage = axes;
                    set(f1,'Name', 'Last Recon');
                end    
                imshow(recent_im, 'Parent', last_montage);
            end
            if length(PS_recon.recon_dat.obj_int) > 1
               for ig = 1:(length(PS_recon.recon_dat.obj_int)-1)
                    intermediate_name = sprintf('%s_int_%.0f_summary.tif',last_saved_file,ig);
                    savesummaryimage3(intermediate_name , ...
                        PS_recon.recon_dat.obj_int{ig},...
                        PS_recon.recon_dat.probe_int{ig},...
                        -rotation_angle);
    %                 SerializeReconInfo(PS_recon.recon_dat));
                    if exec_params.copy_to_monitor_dir 
                        try
                            copyfile(intermediate_name, exec_params.local_monitor_dir);
                        catch
                            fprintf('Something went wrong with copy ... No biggie\n');
                        end
                    end
               end
            end
%             if recipe_var == experi_def.recipe_table.prime_probe_num
%                 prime_recon_fname = last_saved_file;
%             end
    %         %%
        
        catch ME
             warning('\n=== Recipe # %3.2f failed, related to error in === \nFile %s Line: %i \n  id:%s msg:%s \n', ...
                        recipe_var, ME.stack(1).name,  ME.stack(1).line, ME.identifier, ME.message );
        end
        %%
        
        end
    
%             end
    catch pfe
        getReport(pfe,'extended' )
%         error('The error is here (1)');
    end
    
    

    
end    

end_time = now;
exec_time = (end_time - start_time)*24*60*60;
fprintf('Execution Time %.0f\n', exec_time);


catch le
    getReport(le,'extended' )
end

end
%{
General Scheme:
1 - Re-run existing recipe
    -- Calculate centre to inchorerent pedestal ratio used.
    -- See how this compares to that used
    -- Store a recon with the recipe stored with it!
    -- Also store a recon in (complex) MRC format.
2 - Step down pedestal to 2, and use probe and object guess from step 1, 
3 - Step down pedestal to 0, and use probe and object guess from step 2
4 - Step down pedestal to 2, use probe guess from 1 and blank object
5 - Step down pedestal to 0, use probe guess from 4, and object from step 4
---
Investigate DP probe centering and position correction
..

%}

