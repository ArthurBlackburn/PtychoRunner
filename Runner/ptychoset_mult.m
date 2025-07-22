classdef ptychoset_mult < handle
% My main class For collection, storage and processing of ptyhcographic data.
% Contains hangovers from earlier incarnations of data collection and processing workflow, so there
% is likely some dead code in here, though I have tried tu purge it...
%
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
    
        
    properties (SetAccess = public)
        cal_0
        im_posns
        dps
        
        fileObj
        
        current_set
        
        hht
        cam
        
        save_dir
        data_setname
        current_save_path
        
        coms
        recent_received
        
        camera_setup
        BHP
        
        recon_dat
end
        
    methods
       function ps = ptychoset_mult()
       end
       
                    
       function check_and_set_type(ps, set_type)
            allowed_types = {'cal_0','im_posns','dps'};
            if sum(strcmpi(set_type, allowed_types)) == 1
                ps.current_set = set_type;
            else
                error('Invalid data set specified');
            end
        end
        
       function add_data(ps, basename_or_struct, num_range, set_type, append)
                if nargin < 5
                    append = 0;
                end

                if ischar(basename_or_struct)
                    if nargin == 3
                        set_type = 'dps';                        
                    end
                    if nargin == 2
                        num_range = 1;
                    end
                    ps.check_and_set_type(set_type);
                    if isempty(num_range)
                        num_files = 1;
                    else
                        num_files = length(num_range);
                    end
                    filenames = cell(1,num_files);
                    for ii=1:num_files
                        if ~iscell(num_range)
                            if ~isempty(num_range)
                                im_ind = num_range(ii);
                                num_part = num2str(im_ind,'%04.0f');
                            else
                                num_part = [];
                            end
                        else
                            num_part = num_range{ii};
                        end
                        filenames{1,ii} = [basename_or_struct,num_part,'.dm3'];
                    end
                    
                    ps.(ps.current_set).opened_from = basename_or_struct;
                    if append == 0
                        ps.(ps.current_set).suffixes = num_range;
                        
                        [ps.(ps.current_set).ims, ps.(ps.current_set).tags] = read_dm3stack(filenames);
                       
                    else
                        if ~isfield(ps.(ps.current_set),'ims')
                            ps.(ps.current_set).ims = cell(0);
                            ps.(ps.current_set).tags = cell(0);
                            ps.(ps.current_set).suffixes = [];
                        end
                        
                        ps.(ps.current_set).suffixes = [ps.(ps.current_set).suffixes, num_range];
                        last_im = length(ps.(ps.current_set).ims);
                        % grow the cells
                        ps.(ps.current_set).ims = [ps.(ps.current_set).ims, cell(1,num_files)]; 
                        ps.(ps.current_set).tags = [ps.(ps.current_set).tags, cell(1,num_files)];
                        
                        [ps.(ps.current_set).ims((last_im+1):end), ...
                         ps.(ps.current_set).tags((last_im+1):end)] = read_dm3stack(filenames);

                    end
                        % assume that Scale and Units are the same as the
                        % first sets loaded.
                     [ps.(ps.current_set).units_per_pix, ~] = ...
                        ps.(ps.current_set).tags{1}.TheData.Calibrations.Dimensions.Scale;
                        [ps.(ps.current_set).units, ~] = ...
                        ps.(ps.current_set).tags{1}.TheData.Calibrations.Dimensions.Units;
                
                end
                
                if isstruct(basename_or_struct)
                    field_names = fieldnames(basename_or_struct);
                    for f_ind = 1:length(field_names)
                        ps.(field_names{f_ind}) = basename_or_struct.(field_names{f_ind});
                    end
                end
           
       end
              
       function delete(ps)
           if isa(ps.coms, 'instrument')
               fclose(ps.coms);
               delete(ps.coms);
           end
           % Once upon a time this class had functions to control the TEM
           % through Maestro, and used the piBox. 
           % Some commands were passed to Maestro serially and the piBox object was
           % destroyed here. See the archives!        
       end
       
       
       function command_parse(ps)
           cmd_start_ind = strfind(ps.recent_received,'cmd');
           if ~isempty( cmd_start_ind )
               cmd_chain = strsplit(ps.recent_received,':');
%                last_cell = find(strcmpi('end',cmd_chain),1,'first')-1;
               first_cell = find(strcmpi('cmd',cmd_chain),1,'first')+1;
               % allowed commands:
               switch cmd_chain{first_cell}
                   case 'acquire_cal_0'
                       max_shift = str2double(cmd_chain{first_cell + 1});
                       ps.Acquire_Cal_0(max_shift);
                   case 'acquire_fine_pos_reg'
                       FOV = str2double(cmd_chain{first_cell + 1});
                       number_intervals = str2double(cmd_chain{first_cell + 2});
                       n_or_replace = str2double(cmd_chain{first_cell + 3});
                       ps.Acquire_Fine_Pos_Reg(FOV, number_intervals, n_or_replace);
                   case 'acquire_dp_set'
                       grid_step = str2double(cmd_chain{first_cell + 1});
                       num_pts = str2double(cmd_chain{first_cell + 2});
                       rand_frac = str2double(cmd_chain{first_cell + 3});
                       new_set = str2double(cmd_chain{first_cell + 4});
                       ps.Acquire_DP_Set(grid_step, num_pts, rand_frac, new_set);
                   case 'load_dp'
                       set_type = 'dps';
                       newfilename = cmd_chain{first_cell + 1};
                       position = str2num(cmd_chain{first_cell + 2}); %#ok<ST2NM>
                       seq_number = str2double(cmd_chain{first_cell + 3});                       
                       ps.add_to_loadlist(set_type ,newfilename, position', seq_number)
                   case 'test'
                       fprintf('%s\n',cmd_chain{:});
                   otherwise
                       disp('Unknown command received');
               end
           end
       end
       
       function ptycho_ser_rec (ps, obj, event)
           EventType = event.Type;
           name = get(obj, 'Name');
           EventDataTime = event.Data.AbsTime;
           fprintf(['maestro_ser_rec saw', EventType ,'occurring at ' datestr(EventDataTime,13),...
	       ' for the object: ' name '.\n']);
           ps.recent_received = fscanf(obj);           
       end
       
       function setnames(ps, save_dir, data_setname)
           ps.save_dir = save_dir;
           ps.data_setname = data_setname;
           ps.current_save_path = fullfile(ps.save_dir, ps.data_setname );
       end
       
       function okay = MakeDirs( ps )
        % Makes directories for storing ptycho data if it they are not already made

        dirs_to_make = {...
            'cal_0',...
            'camera_setup',...
            'dps',...
            'im_posns',...
            'recons',...
            'recon_dat'};

        okay = cell(length(dirs_to_make),2);

        for ci = 1:length(dirs_to_make)
            if ~exist(fullfile(ps.save_dir, ps.data_setname, dirs_to_make{ci} ), 'dir')
                [okay{ci,1}, okay{ci,2}] = ...
                mkdir(fullfile(ps.save_dir, ps.data_setname), dirs_to_make{ci});                
            end
        end
        
       end
       
       function Acquire_Cal_0(ps, max_shift)
           % setup the camera, make sure exposure is okay in Maestro first:
           savepath = [fullfile(ps.current_save_path,'cal_0'),'\'];
           filename_base = 'cal_0_';
           ps.cam.Settings.Acquire.Save = 'on';
           ps.cam.Settings.Acquire.Display = 'new';
           ps.cam.Settings.Acquire.SavePath = savepath;
           if ~isprop(ps,'BHP') || ~isfield(ps.BHP,'cal_cent_point') || isempty(ps.BHP.cal_cent_point)     
               ps.BHP.cal_cent_point = ps.hht.Coils.BH;
           end
           start_position = ps.BHP.cal_cent_point;
           ps.cam.Upper.ImageName = [filename_base,'0_0'];
           ps.cam.Upper.AcquireImage;
           ps.hht.Coils.BH = start_position + [max_shift, 0];
           ps.cam.Upper.ImageName = [filename_base,'1_0'];
           ps.cam.Upper.AcquireImage;
           ps.hht.Coils.BH = start_position + [00, max_shift];
           ps.cam.Upper.ImageName =  [filename_base,'0_1'];
           ps.cam.Upper.AcquireImage;
           ps.hht.Coils.BH = start_position; 
                  
       end
       
       function new_or_replace_existing_calset(ps)
               savepath = [fullfile(ps.current_save_path,'cal_0'),'\'];
               filename_base = 'cal_0_';
               ps.cal_0.opened_from = [savepath, filename_base];
               ps.cal_0.suffixes = {'0_0','1_0','0_1'};
               ps.cal_0.ims = cell(0);
               ps.analyze_cal_0;
       end
       
       function create_cfits_on_deltas(ps, x_act, y_act, del_x, del_y)
            ft = fittype( 'poly11' );
            [xData, yData, zData] = prepareSurfaceData( x_act, y_act, del_x );
            act_XY_to_DELX = fit( [xData, yData], zData, ft );
            [Pact_XY_to_DELX, a] = ps.conv_sfits_to_poly(act_XY_to_DELX);
%             plot( act_XY_to_DELX, [xData, yData], zData );
            [xData, yData, zData] = prepareSurfaceData( x_act, y_act, del_y );
            act_XY_to_DELY = fit( [xData, yData], zData, ft );
            [Pact_XY_to_DELY, b] = ps.conv_sfits_to_poly(act_XY_to_DELY);
%             plot( act_XY_to_DELY, [xData, yData], zData )
            
            ps.im_posns.act_XY_to_deltaXY = @(cd) [Pact_XY_to_DELX(cd(1), cd(2));...
                                                   Pact_XY_to_DELY(cd(1), cd(2))];
            ps.im_posns.act_XY_to_deltaXY_coefs = {a, b};
            
            x_req = x_act - del_x;
            y_req = y_act - del_y;
            
            [xData, yData, zData] = prepareSurfaceData( x_req, y_req, del_x );
            req_XY_to_DELX = fit( [xData, yData], zData, ft );
            [Preq_XY_to_DELX, c] = ps.conv_sfits_to_poly(req_XY_to_DELX);
%             plot( req_XY_to_DELX, [xData, yData], zData );


            [xData, yData, zData] = prepareSurfaceData( x_req, y_req, del_y );
            req_XY_to_DELY = fit( [xData, yData], zData, ft );
            [Preq_XY_to_DELY, d] = ps.conv_sfits_to_poly(req_XY_to_DELY);
%             plot( req_XY_to_DELY, [xData, yData], zData );
            ps.im_posns.req_XY_to_deltaXY = @(cd) [Preq_XY_to_DELX(cd(1), cd(2));...
                                                   Preq_XY_to_DELY(cd(1), cd(2))];
            ps.im_posns.req_XY_to_deltaXY_coefs = {c, d};
       end
       
       function inputs = conv_ptychoset_lsmql(ps)
            if isa(ps.dps.ims_proc, 'uint16')
                inputs.diffraction = single(ps.dps.ims_proc);
            else
                inputs.diffraction = ps.dps.ims_proc;
            end
            if isfield(ps.dps,'masks')
                inputs.masks = ~ps.dps.masks;
            end
            if isfield(ps.dps,'edgezeros')
                inputs.edgezeros =  ps.dps.edgezeros;
            end
            
            inputs.object_orig = ps.recon_dat.object_0.'; % testing to see if this resolves the subregion problem with using Matlab ePIE implementation.
            inputs.probe_orig = ps.recon_dat.probe_pf.'; % testing to see if this resolves the subregion problem with using Matlab ePIE implementation.
            inputs.probe_positions_0 = [ps.dps.dpPosR',ps.dps.dpPosC'];
            inputs.probe_positions_0 = inputs.probe_positions_0  - mean(inputs.probe_positions_0);
            inputs.reconstruct_ind{1} = 1:length(inputs.probe_positions_0);
            inputs.Np_p = size(inputs.probe_orig);
            inputs.Np_o = size(inputs.object_orig);
            inputs.Npos = length(inputs.probe_positions_0);
            inputs.object{1} = inputs.object_orig;
            inputs.probe{1} = inputs.probe_orig;
            inputs.probe_positions = [];
            inputs.probe_support = [];
            inputs.modes = [];
            inputs.probe_evolution = [];
       end
       
 
       function recreate_im_posn_cfit(ps)
           a = ps.create_sfits_from_coefs(ps.im_posns.act_XY_to_deltaXY_coefs{1});
           b = ps.create_sfits_from_coefs(ps.im_posns.act_XY_to_deltaXY_coefs{2});
           c = ps.create_sfits_from_coefs(ps.im_posns.req_XY_to_deltaXY_coefs{1});
           d = ps.create_sfits_from_coefs(ps.im_posns.req_XY_to_deltaXY_coefs{2});
           
           ps.im_posns.act_XY_to_deltaXY = @(cd) ...
               [a(cd(1), cd(2)) ;...
                b(cd(1), cd(2)) ];
            
           ps.im_posns.req_XY_to_deltaXY = @(cd) ...
               [c(cd(1), cd(2)) ;...
                d(cd(1), cd(2)) ]; 
       end

       function posn_fit_method(ps)
           if ~isfield(ps.im_posns,'use_cfit') || (ps.im_posns.use_cfit == 0)
               ps.im_posns.posn_method = @(go_to_posn, base_BH) ps.BH_from_xy_cal_0(go_to_posn, base_BH);
           else
           % Use cfit
               ps.im_posns.posn_method = @(go_to_posn, base_BH)...
                   ps.BH_from_xy_cal_0( ps.request_from_actual(go_to_posn) , base_BH);
           end
       end
       
       function use_im_posn_cfit(ps, on_or_off)
           switch on_or_off
               case 'on'
                   ps.im_posns.use_cfit = 1;
               case 'off'
                   ps.im_posns.use_cfit = 0;
               otherwise
                   disp('invalid option');
           end
       end
       
       function req_position = request_from_actual(ps, act_pos)
           req_position = act_pos - ps.im_posns.act_XY_to_deltaXY(act_pos);
           % as delta = act_pos - requested pos
           % thus req_pos = act_ps - delta;
       end
       
       function act_position = actual_from_request(ps, req_pos)
           act_position = req_pos + ps.im_posns.req_XY_to_deltaXY(req_pos);
       end
       
       function Acquire_Fine_Pos_Reg(ps, FOV, number_intervals, new_or_replace)
           setname = 'im_posns';
           ps.posn_fit_method();
           %%%%----
           posn_method = ps.im_posns.posn_method;
           %%%%----
           rand_scaler = 0;
           ps.Acquire_Fine_Grid_of_Ims(FOV, number_intervals, rand_scaler, setname, posn_method, new_or_replace);
       end
       
       function Acquire_DP_Set(ps, grid_step, num_pts, rand_fraction, new_or_replace, remote_load)
           if nargin < 6
               remote_load = 0;
           end
           
           FOV = num_pts*grid_step;
           number_intervals = num_pts - 1;
           rand_scaler = rand_fraction * grid_step;
           setname = 'dps';
           ps.posn_fit_method();
           % can later change this to the function fit from previous stage
           
           posn_method = ps.im_posns.posn_method;
           %%%--
           ps.Acquire_Fine_Grid_of_Ims(FOV, number_intervals, rand_scaler, setname, ...
               posn_method, new_or_replace, remote_load)
           ps.hht.Coils.BH = ps.BHP.ptycho_safe_place;
       end
       
       function shift_dps_by_CofM_of_first(ps, show_ops)
           cent_pos_1 = FindCentre(ps.dps.ims{1});
           ps.dps.ims_proc = cell(size(ps.dps.ims));
           for n = 1:length(ps.dps.ims)
                ps.dps.ims_proc{n} = abs(Recentre(ps.dps.ims{n},cent_pos_1));
                if show_ops
                    imagesc(ps.dps.ims{n});
                    drawnow;
                end
           end
       end
       
       function crop_dps(ps, operate_on, cropped_size)
           base_size = size(ps.dps.(operate_on){1});
           r_rng = ((base_size(1)-cropped_size)/2 + 1):...
                   ((base_size(1)-cropped_size)/2 + cropped_size) ;
           c_rng = ((base_size(2)-cropped_size)/2 + 1):...
                   ((base_size(2)-cropped_size)/2 + cropped_size) ;    
           for n = 1:length(ps.dps.(operate_on))
                ps.dps.ims_proc{n} = ps.dps.(operate_on){n}(r_rng, c_rng);
           end
       end
       
       function get_dp_positions_fromtag(ps, pos_det_method)
            % reference things to this first position
            BH_1 = ...
               [ps.dps.tags{1}.Microscope.Coils.BH.x;...
                ps.dps.tags{1}.Microscope.Coils.BH.y];
            if nargin < 2
                pos_det_method = 'simple';
            end
            switch pos_det_method
                case 'simple'
                    pos_func = @(BH_dp) ps.xy_from_BH_cal_0(BH_dp, BH_1); 
                case 'cift'
                    pos_func = @(BH_dp) ps.actual_from_request(ps.xy_from_BH_cal_0(BH_dp, BH_1));
                otherwise
                    error('unknown method');
            end
            
            ps.dps.scan_posns_fromtags = zeros(2,length(ps.dps.tags));

                for ii = 2:length(ps.dps.tags)
                   BH_dp = ...
                       [ps.dps.tags{ii}.Microscope.Coils.BH.x;...
                        ps.dps.tags{ii}.Microscope.Coils.BH.y];                        
                   ps.dps.scan_posns(:,ii) = ...
                          pos_func(BH_dp); 
                end
       end
       
       function get_large_area_zoom_image(ps, imname)
           pos_start = 'image_safe_place';
           pos_exp = 'image_exposure_pos';
           pos_end = 'image_safe_place';
           C1_start = 'C1_image';
           C1_end = [];
           apt_start = 'Apt1';
           apt_end = 'Apt1';
           camera_set = 'zoom_large';
           if ~isempty(imname)
               save_acquisition = 'on';
           else
               save_acquisition = 'off';
           end
           ps.get_image(imname, ...
               pos_start, pos_exp, pos_end,...
               C1_start, C1_end,...
               apt_start, apt_end, ...
               camera_set, save_acquisition)
       end
       
       function get_small_area_zoom_image(ps, imname)
           pos_start = 'ptycho_safe_place';
           pos_exp = 'ptycho_safe_place';
           pos_end = 'ptycho_safe_place';
           C1_start = 'C1_ptycho';
           C1_end = [];
           apt_start = 'Apt3'; %% need to make variable!
           apt_end = 'Apt3'; %% need to make variable!
           camera_set = 'zoom_small';
           if ~isempty(imname)
               save_acquisition = 'on';
           else
               save_acquisition = 'off';
           end
           ps.get_image(imname, ...
               pos_start, pos_exp, pos_end,...
               C1_start, C1_end,...
               apt_start, apt_end, ...
               camera_set, save_acquisition)
       end

       function get_diff_image(ps, imname)
           pos_start = 'ptycho_safe_place';
           pos_exp = 'ptycho_safe_place';
           pos_end = 'ptycho_safe_place';
           C1_start = 'C1_ptycho';
           C1_end = [];
           apt_start = 'Apt3'; %% need to make variable!
           apt_end = 'Apt3'; %% need to make variable!
           camera_set = 'diff';
           if ~isempty(imname)
               save_acquisition = 'on';
           else
               save_acquisition = 'off';
           end
           ps.get_image(imname, ...
               pos_start, pos_exp, pos_end,...
               C1_start, C1_end,...
               apt_start, apt_end, ...
               camera_set, save_acquisition)
       end
       
       
       function get_image(ps, imname, ...
               pos_start, pos_exp, pos_end,...
               C1_start, C1_end,...
               apt_start, apt_end, ...
               camera_set, save_acquisition)
           
           if nargin < 10
               save_acquisition = 'on';
           end
           
           if ischar(pos_start)
               pos_start = ps.BHP.(pos_start);
           end
           
           if ischar(pos_exp)
               pos_exp = ps.BHP.(pos_exp);
           end
           
           if ischar(pos_end)
               pos_end = ps.BHP.(pos_end);
           end
           
           savepath = [fullfile(ps.current_save_path,'camera_setup'),'\'];
           ps.cam.Settings.Acquire.Save = save_acquisition;
           ps.cam.Settings.Acquire.Display = 'new';
           ps.cam.Settings.Acquire.SavePath = savepath;
           
           ps.hht.Coils.BH = pos_start;
           ps.hht.Lenses.C1 = ps.BHP.(C1_start);
           ps.hht.Apertures.COND = apt_start;
           ps.set_for_imaging(camera_set);
           if ~isempty(imname)
               ps.cam.Upper.ImageName = imname;
           else
               ps.cam.Upper.ImageName = 'Temp';
           end
           ps.hht.Coils.BH = pos_exp;
           ps.cam.Upper.AcquireImage();
           ps.hht.Coils.BH = pos_end;
           if ~isempty(C1_end)
               ps.hht.Lenses.C1 = ps.BHP.(C1_end);
           end
           ps.hht.Apertures.COND = apt_end;        
           
       end
       
       function set_for_imaging(ps, imode_in ,image_name, imexposure, binning, imfrac)
           switch imode_in
               case 'zoom_large'
                   imode = 'zoom_large';
               case 'zoom_small'
                   imode = 'zoom_small';
               case 'diff'
                   imode = 'diff';
               otherwise
                   error('mode not recognized')
           end
           
           if nargin < 6
               imfrac = ps.camera_setup.(imode).imfrac;
           else
               ps.camera_setup.(imode).imfrac = imfrac;
           end
           ps.set_camera_portion(imfrac);
           if nargin < 5
               binning = ps.camera_setup.(imode).binning;
           else
               ps.camera_setup.(imode).binning = binning;
           end
           ps.cam.Upper.Binning.x = binning ;
           ps.cam.Upper.Binning.y = binning ;
           if nargin < 4
               imexposure = ps.camera_setup.(imode).imexposure;
           else
               ps.camera_setup.(imode).imexposure = imexposure;
           end
           if nargin < 3
               image_name = [];
           end
           ps.cam.Upper.Exposure = imexposure;
           ps.camera_setup.last_mode_set = imode;
           if ~isempty(image_name)
                ps.cam.Upper.ImageName = image_name;
                ps.cam.Upper.AcquireImage();  
           end           
       end
       
       function set_camera_portion(ps, cam_frac)
            cam_sects = [...
                0, 0, 2048, 2048;
                512,512,1536,1536;
                768, 768, 1280, 1280];
            top = cam_sects(:,1);
            left = cam_sects(:,2);
            bottom = cam_sects(:,3);
            right = cam_sects(:,4);
            cam_tbl = table(top, left, bottom, right, 'RowNames', {'full','half','quarter'});
            ps.cam.Upper.Area.top = cam_tbl{cam_frac, 'top'};
            ps.cam.Upper.Area.left = cam_tbl{cam_frac, 'left'};
            ps.cam.Upper.Area.bottom = cam_tbl{cam_frac, 'bottom'};
            ps.cam.Upper.Area.right = cam_tbl{cam_frac, 'right'};             
       end
       
       function Acquire_Fine_Grid_of_Ims(ps, FOV, number_intervals, rand_scaler, ...
               setname, posn_method, new_or_replace, send_remote_load_req)
           
           if nargin < 8
               send_remote_load_req = 0;
           end
           
           savepath = [fullfile(ps.current_save_path,setname),'\'];
           filename_base = [setname,'_'];
           
           ps.cam.Settings.Acquire.Save = 'on';
           ps.cam.Settings.Acquire.Display = 'off';
           ps.cam.Settings.Acquire.SavePath = savepath;
           
           if ~isprop(ps,'BHP') || ~isfield(ps.BHP,'ptycho_start_point') || isempty(ps.BHP.ptycho_start_point)     
               ps.BHP.ptycho_start_point = ps.hht.Coils.BH;
           end
           
           base_BH = [ps.BHP.ptycho_start_point.x;...
                      ps.BHP.ptycho_start_point.y]; 
                  
           pattern_type = 'meander_x_scan';
           step_size = FOV/number_intervals;
           [xo, yo] = meshgrid(0:step_size:FOV);
           scan_seq = ps.sequence_points(xo, pattern_type);
           scan_posns = zeros(2,length(scan_seq));
                      
           if rand_scaler ~= 0
               % originally this was as below - but noted with MH the
               % asymmetry. Thus changed Jan 2017:
%                rand_x = rand(size(xo));
%                xo = xo + rand_scaler* rand_x;
               rand_x = rand(size(xo)) - 0.5;
               rand_y = rand(size(xo)) - 0.5;
               xo = xo + 2*(rand_scaler* rand_x);
               yo = yo + 2*(rand_scaler* rand_y);
           end
           
           req_postions = zeros(numel(xo),2);
           
           suffixes = 1:1:numel(xo);
           for im_number = suffixes
               go_to_posn = [xo(scan_seq(im_number));...
                   yo(scan_seq(im_number))];
               req_postions(im_number,:) = go_to_posn';
               ps.(setname).opened_from.req_postions = req_postions;
               scan_posns(:,im_number) = go_to_posn;
               go_to_BH = posn_method(go_to_posn, base_BH);
               ps.hht.Coils.BH = go_to_BH';
               sprintf('Now at x: %6.3f, y %6.3f',go_to_posn);
               num_part = num2str(im_number ,'%04.0f');
               imname = [filename_base,num_part];
               ps.cam.Upper.ImageName = imname;
               ps.cam.Upper.AcquireImage;
               disp(['Saved to file ',imname]);
               if send_remote_load_req
                   fprintf(ps.coms,['cmd:load_dp:',imname,':',...
                                      num2str(go_to_posn'),':'...
                                       num2str(im_number),':end']);
               end
           end
           save([savepath,filename_base,'req_posns.mat'] ,'req_postions');
           % put beam back to home
           if new_or_replace
               ps.(setname).opened_from = [savepath, filename_base];
               ps.(setname).suffixes = suffixes;
               ps.(setname).scan_posns = scan_posns;
               ps.(setname).base_BH = base_BH;
               ps.(setname).pattern_type = pattern_type;
               ps.(setname).ims = cell(0);
           end
       end
       
       function Analyze_Fine_Pos_Reg(ps, method_str)
               % load data if needed:
               if isempty(ps.im_posns.ims)
                   ps.add_data(ps.im_posns.opened_from, ps.im_posns.suffixes,'im_posns');
               end
               if nargin < 2 || isempty(method_str)
                   method_str = 'xcorr';
               end
               ps.current_set = 'im_posns';
               reference_index = 1;
               ps.find_offsets(method_str, reference_index);
               ps.sort_offsets;
%                ps.im_posns.sorted_offsets = ps.im_posns.raw_offsets;
               ps.scale_offsets;
               ps.im_posns.calced_offsets = ...
               [ps.im_posns.scaled_offsets.x;ps.im_posns.scaled_offsets.y];
               if ~isfield( ps.im_posns,'scan_posns')
                   numpts = length(ps.im_posns.suffixes);
                   ps.im_posns.scan_posns = zeros(2, numpts);
                   BH_1 = ...
                           [ps.im_posns.tags{1}.Microscope.Coils.BH.x;...
                            ps.im_posns.tags{1}.Microscope.Coils.BH.y];
                   for ii = 2:numpts
                       BH_im = ...
                           [ps.im_posns.tags{ii}.Microscope.Coils.BH.x;...
                            ps.im_posns.tags{ii}.Microscope.Coils.BH.y];                        
                       ps.im_posns.scan_posns(:,ii) = ...
                           ps.xy_from_BH_cal_0(BH_im, BH_1);                             
                   end
               end
               % Delta = Measured - Requested, so to get to desired point more accurately
               % Request should be (Desired Point - Delta), where Delta is
               % found from looking at map of Delta vs measured position.
               % similarly if aiming to determine poistion from coil
               % values, accuarte poistion comes from looking at simple
               % requested value + delta, where delta comes from map with
               % simple requested value
               ps.im_posns.offset_deltas = ...
                   ps.im_posns.calced_offsets -  ...
                   ps.im_posns.scan_posns;           
       end
       
       function plot_im_posn_ReqToAct(ps)
            figure; 
            plot(ps.im_posns.scan_posns(1,:),...
                 ps.im_posns.scan_posns(2,:),'r');
            hold on
            plot(ps.im_posns.calced_offsets(1,:),...
                 ps.im_posns.calced_offsets(2,:),'b');
            axis equal
            xlabel('X (nm)'); ylabel('Y (nm)');
            legend('Requested Input', 'Actual from TEM', 'Location','NorthOutside',...
                'Orientation','horizontal');
            legend boxoff
       end
       
       function [dp_posns, dp_indices] = create_subset(ps, central_RCsize, operate_on)
            if nargin < 3
                operate_on = 'scan_posns';
            end
            
            if nargin < 2
                central_RCsize = [];
            end
            
            xpts = ps.dps.(operate_on)(1,:);
            ypts = ps.dps.(operate_on)(2,:);
            
            dp_indices = 1:size(ps.dps.(operate_on),2);
            
            if ~isempty(central_RCsize)
                % Then use this sized region from the central pattern:
                % Check that the selected size is smaller than the source:
                if ~(isfield(ps.dps,'attrs') && isfield(ps.dps.attrs,'meshParams') && isfield(ps.dps.attrs.meshParams,'shape'))
                    error('Mesh of positions has unknown shape.\n');
                end     
                if any(central_RCsize > ps.dps.attrs.meshParams.shape)
                    error('Entered subregion Size is greater than the mesh size.\n');
%                              'Reverting to manual selection of points from GUI'];
%                   central_RCsize = []; % replace the above 'error' by
%                   'warning' if so desired. Doing so would not be good for
%                   batch processing though.
                end           
                rownums = 1:ps.dps.attrs.meshParams.shape(1);
                colnums = 1:ps.dps.attrs.meshParams.shape(2);
                [C, R] = meshgrid(colnums, rownums);
                r1 = round((ps.dps.attrs.meshParams.shape(1)/2) - (central_RCsize(1)/2));
                r2 = r1 + (central_RCsize(1)-1);
                c1 = round((ps.dps.attrs.meshParams.shape(2)/2) - (central_RCsize(2)/2));
                c2 = c1 + (central_RCsize(2)-1);
                inds_in_bnds = ( r1 <= R) & (R <= r2) & ...
                               ( c1 <= C) & (C <= c2);
                % if the above is read off using (:) operator it goes down
                % the columns. If we are use this to access the order in
                % normal raster scan, we want it to be read along rows,
                % fromt the top. i.e. need to transpose:
                inds_in_bnds = inds_in_bnds.';
                
            end
            
            if nargin == 1 || isempty(central_RCsize)
                pts_sel = figure;
                plot(xpts, ypts, 'x'); axis equal;
                disp('Select Sub Region');
                rect_region = getrect(pts_sel);
%                 close(pts_sel);
                inds_in_bnds = ...
                    (rect_region(1) < xpts) & ( xpts < (rect_region(1) + rect_region(3))) & ...
                    (rect_region(2) < ypts) & ( ypts < (rect_region(2) + rect_region(4)));
            end
         
            dp_posns = ps.dps.(operate_on)(:,inds_in_bnds);
            dp_indices = dp_indices(inds_in_bnds);
            
%             x_bnd = rect_region(1);
%             y_bnd = rect_region(2);
%             offsets.x = raw_offsets.x;
%             offsets.y = raw_offsets.y;
%             offsets.x(offsets.x >= x_bnd) = -(dim2 - offsets.x(offsets.x >= x_bnd));
%             offsets.y(offsets.y >= y_bnd) = -(dim1 - offsets.y(offsets.y >= y_bnd));
%             offsets.x = offsets.x - offsets.x(1);
%             offsets.y = offsets.y - offsets.y(1);
           
       end
       
       function [x_actual, y_actual, del_x, del_y] = plot_deltas(ps)
            x_actual = ps.im_posns.calced_offsets(1,:);
            y_actual = ps.im_posns.calced_offsets(2,:);
%             x = ps.im_posns.scan_posns(1,:);
%             y = ps.im_posns.scan_posns(2,:);
            del_x = ps.im_posns.offset_deltas(1,:);
            del_y = ps.im_posns.offset_deltas(2,:);
            figure; 
            plot3(x_actual, ...
                  y_actual, ...
                  ps.im_posns.offset_deltas(1,:),'r'); 
            xlabel('X actual (nm)'); ylabel('Y actual(nm)'); zlabel('Detla X / Y (nm)');
            hold on;
            plot3(x_actual, ...
                  y_actual, ...
                  ps.im_posns.offset_deltas(2,:),'b'); 
            legend('Delta X (nm)', 'Delta Y (nm)', 'Location','NorthOutside',...
                'Orientation','horizontal');
            legend boxoff
        end
       
       function ListenForNewFiles(ps, pathToWatch)
           ps.fileObj = System.IO.FileSystemWatcher(pathToWatch);
           ps.fileObj.EnableRaisingEvents = true;
           addlistener(ps.fileObj, 'Changed', @ps.NewFileDet);
       end
       
       function StopListenForNewFiles(ps)
           ps.fileObj.EnableRaisingEvents = false;
           
       end
       
       function NewFileDet(ps,~,evt)
           % listener action if comboBox changes
           new_file = char(evt.FullPath.ToString());
           disp(['New file Detected: ',new_file]);
           ps.check_and_load(ps.current_set, new_file);
       end
                     
       function LoadFileFromList_on_Arrival(ps, set_type)
%            ps.ListenForNewFiles(pathToWatch, @ps.check_and_load);
           pathToWatch = fullfile(ps.current_save_path,set_type);
           ps.ListenForNewFiles(pathToWatch);
           ps.current_set = set_type;
           
       end   
       
%        function LoadFileFromList_AtInterval(ps, est_type, interval)
%            
%        end
%        
       function check_and_load(ps, set_type, new_file)
           pathToWatch = fullfile(ps.current_save_path,set_type);
           
           if isfield(ps.(set_type), 'filelist_to_load')
               num_that_should_be_loaded = ...
                   ps.(set_type).filelist_to_load.num_files_in_list;
           else
               num_that_should_be_loaded = [];
           end
           
           if  ~isempty(num_that_should_be_loaded) && ...
               (num_that_should_be_loaded ~= sum(ps.(set_type).filelist_to_load.loaded))
               lookfor = find ( ~ps.(set_type).filelist_to_load.loaded(1:num_that_should_be_loaded) );
           else
               disp('New file appeared that is not on the load list');
               lookfor = [];
           end
           if ~isempty(new_file)
           [~,NAME,EXT] = fileparts(new_file);
           else
               NAME = [];
           end
           for li = 1:length(lookfor)
               list_ind = lookfor(li);
               % is it the file just flagged up?
               looking_for_file = ps.(set_type).filelist_to_load.filenames{list_ind};
               if strcmpi(looking_for_file, NAME) || ... % or is perhaps already in this directory, perhaps missed before
                  exist(fullfile(pathToWatch, [ps.(set_type).filelist_to_load.filenames{list_ind},'.dm3']),'file')
                   
                   ps.add_data( [pathToWatch,'\'], ...
                          ps.(set_type).filelist_to_load.filenames(list_ind), set_type, 1 );
                      
                   disp(['Loaded File ', ps.(set_type).filelist_to_load.filenames{list_ind},...
                       ' into ',set_type]);
                   
                   ps.(set_type).filelist_to_load.loaded(list_ind ) = 1;
               end
           end
       end
                     
       function analyze_cal_0(ps, method_str)
           % load data if needed:
           if isempty(ps.cal_0.ims)
               ps.add_data(ps.cal_0.opened_from, ps.cal_0.suffixes,'cal_0');
           end
           if nargin < 2 || isempty(method_str)
               method_str = 'xcorr';
           end
           ps.current_set = 'cal_0';
           % sort out which image is which:
           if ischar(ps.(ps.current_set).suffixes)
               BH_base_ind = find(strcmpi('0_0',ps.(ps.current_set).suffixes),1);
               BHx_plus_image_ind = find(strcmpi('1_0',ps.(ps.current_set).suffixes),1);
               BHy_plus_image_ind = find(strcmpi('0_1',ps.(ps.current_set).suffixes),1);
           else
               BH_base_ind = 1;
               BHx_plus_image_ind = 2;
               BHy_plus_image_ind = 3;
           end
           ps.find_offsets(method_str, BH_base_ind);
           BHx_plus_image_shift = [...
               ps.(ps.current_set).raw_offsets.x(BHx_plus_image_ind),...
               ps.(ps.current_set).raw_offsets.y(BHx_plus_image_ind)];
           BHy_plus_image_shift = [...
               ps.(ps.current_set).raw_offsets.x(BHy_plus_image_ind),...
               ps.(ps.current_set).raw_offsets.y(BHy_plus_image_ind)];  
           
           meas.x_theta = atan2(BHx_plus_image_shift(2), BHx_plus_image_shift(1));
           meas.y_theta = atan2(BHy_plus_image_shift(2), BHy_plus_image_shift(1));
           [meas.inc_theta, meas.flip_y] = IncAngle(meas.x_theta, meas.y_theta, 0);
           meas.BH_base_point = [...
               ps.(ps.current_set).tags{BH_base_ind}.Microscope.Coils.BH.x,...
               ps.(ps.current_set).tags{BH_base_ind}.Microscope.Coils.BH.y];
           meas.BHx_current_shift = ...
               ps.(ps.current_set).tags{BHx_plus_image_ind}.Microscope.Coils.BH.x - ...
               meas.BH_base_point(1);
           meas.BHy_current_shift = ...
               ps.(ps.current_set).tags{BHy_plus_image_ind}.Microscope.Coils.BH.y - ...
               meas.BH_base_point(2);
           meas.BHx_measured_shift = ps.(ps.current_set).units_per_pix * ...
               sqrt(sum(BHx_plus_image_shift.^2));
           meas.BHy_measured_shift = ps.(ps.current_set).units_per_pix * ...
               sqrt(sum(BHy_plus_image_shift.^2));
           
           meas.lin_units = ps.(ps.current_set).units;
           meas.ang_units = 'deg';
           ang_scale = 180/pi;
           meas.x_theta = meas.x_theta * ang_scale;
           meas.y_theta = meas.y_theta * ang_scale;
           meas.inc_theta = meas.inc_theta * ang_scale;
           ps.(ps.current_set).meas = meas; 
           [ps.(ps.current_set).meas.xy_out_T_BH_in,...
            ps.(ps.current_set).meas.bh_out_T_xy_in] = ps.get_3pt_transform;
           
       end
       
       function save(ps,whichpart,basename)
           warning('off','MATLAB:structOnObject');
           
           if nargin < 3
               basename = [];
           end
           
           comps = ls_def_parts(ps, whichpart);
           
           if isempty(basename)
               stem_name = ps.(comps.data_set).opened_from;
           else
               stem_name = basename;                       
           end
           
           save_struct = comps.save_struct(ps); %#ok<NASGU>
           s_name = [stem_name,comps.file_suffix];
           save(s_name,'-struct','save_struct');
       end
       
       function comps = ls_def_parts(ps, whichpart)
            switch whichpart
               case 'cal_0_meas'
                   comps.data_set= 'cal_0';
                   comps.load_struct = 'meas';
                   comps.save_struct = @(ps) struct(ps.(comps.data_set).(comps.load_struct));                    
                   comps.file_suffix = '_cal_0_meas.mat';                  
               case 'cal_0_ims'
                   comps.data_set= 'cal_0';
                   comps.load_struct = {'ims','tags'};
                   comps.save_struct = @(ps) ...
                       struct(comps.load_struct{1}, {ps.(comps.data_set).(comps.load_struct{1})},...
                              comps.load_struct{2}, {ps.(comps.data_set).(comps.load_struct{2})});                      
                   
                   comps.file_suffix = '_cal_0_ims.mat';  
               case 'posn_ims'
                   comps.data_set= 'im_posns';
                   comps.load_struct = {'ims','tags'};
                   comps.save_struct = @(ps) ...
                       struct(comps.load_struct{1}, {ps.(comps.data_set).(comps.load_struct{1})},...
                              comps.load_struct{2}, {ps.(comps.data_set).(comps.load_struct{2})});                      
                   
                   comps.file_suffix = '_imposn_ims.mat';  
               case 'posn_ims_ana'
                   comps.data_set= 'im_posns';
                   if isfield(ps.im_posns, 'tags')
                       comps.save_struct = @(ps) rmfield(ps.im_posns,{'ims','tags'}); 
                   else
                       comps.save_struct = @(ps) struct(ps.(comps.data_set));
                   end
                   comps.load_struct = [];
                   comps.file_suffix = '_imposn_ana.mat';
               case 'posn_ims_cfit'
                   comps.data_set= 'im_posns';
                   comps.load_struct = {'act_XY_to_deltaXY_coefs','req_XY_to_deltaXY_coefs'};
                   comps.save_struct = @(ps) ...
                       struct(comps.load_struct{1}, {ps.(comps.data_set).(comps.load_struct{1})},...
                              comps.load_struct{2}, {ps.(comps.data_set).(comps.load_struct{2})}); 
                   comps.file_suffix = '_imposn_cfit.mat';                   
               case 'dp_ims_raw'                   
                   comps.data_set= 'dps';
                   if isfield(ps.dps,'ims_proc')
                       comps.save_struct = @(ps) rmfield(ps.dps, 'ims_proc');
                   else
                       comps.save_struct = @(ps) ps.dps;
                   end
                   comps.load_struct = [];
                   comps.file_suffix = '_dps_ims_raw.mat';
               case 'dp_ims_proc'
                   comps.data_set= 'dps';
                   comps.load_struct = {'ims_proc','scan_posns'}; % need to fix for when scan pos is saved elsewhere...
                   comps.save_struct = @(ps) ...
                       struct(comps.load_struct{1}, {ps.(comps.data_set).(comps.load_struct{1})},...
                              comps.load_struct{2}, {ps.(comps.data_set).(comps.load_struct{2})}); 
                   comps.file_suffix = '_dps_ims_proc.mat';
               case 'dp_scan_posns'
                   comps.data_set= 'dps';
                   comps.load_struct = {'scan_posns'};
                   comps.save_struct = @(ps) ...
                       struct(comps.load_struct{1}, {ps.(comps.data_set).(comps.load_struct{1})});
                   comps.file_suffix = '_dps_scan_posns.mat';
              case 'dp_scan_posns_fromtags'
                   comps.data_set= 'dps';
                   comps.load_struct = {'scan_posns_fromtags'};
                   comps.save_struct = @(ps) ...
                       struct(comps.load_struct{1}, {ps.(comps.data_set).(comps.load_struct{1})});
                   comps.file_suffix = '_dps_scan_posns_fromtsgs.mat';
               case 'recon_results'
                   comps.data_set= 'recon_dat';
                   comps.load_struct = {'probe_rec','object_rec','recon_posns','recon_pixel_size'};
                   comps.save_struct = @(ps) ...
                       struct(comps.load_struct{1}, {ps.(comps.data_set).(comps.load_struct{1})},...
                              comps.load_struct{2}, {ps.(comps.data_set).(comps.load_struct{2})},...
                              comps.load_struct{3}, {ps.(comps.data_set).(comps.load_struct{3})},...
                              comps.load_struct{4}, {ps.(comps.data_set).(comps.load_struct{4})}); 
                   comps.file_suffix = '_recon_dat.mat';
                    
               otherwise
                   error('unknown type');
           end
        end
       
       function load(ps,whichpart,basename)
           if nargin < 3
               basename = [];
           end
           
           comps = ls_def_parts(ps, whichpart);
           
           if isempty(basename)
               stem_name = ps.(comps.data_set).opened_from;
           else
               stem_name = basename;                       
           end
           
           f_name = [stem_name,comps.file_suffix];
           if ~iscell(comps.load_struct)
               if ~isempty(comps.load_struct)
                   ps.(comps.data_set).(comps.load_struct) = load(f_name);
               else
                   % should put something to save removed fields then add
                   % them back in. Add extra field to comps
                   % ('removed_fields' for example).
                   ps.(comps.data_set) = load(f_name);
               end
           else
               loaded = load(f_name);
               for ii = 1:length(comps.load_struct)
                   ps.(comps.data_set).(comps.load_struct{ii}) = ...
                       loaded.(comps.load_struct{ii});
               end
           end
         
       end
       
       function xy = xy_from_BH_cal_0(ps, BH, BH_at_xyorigin)
            if nargin < 3
                BH_at_xyorigin = ps.cal_0.meas.BH_base_point';
            end
            Tmat = ps.cal_0.meas.xy_out_T_BH_in;
            coords_full = Tmat*[BH,BH_at_xyorigin; ones(1,1+size(BH,2))];
            xy1 = coords_full(1:2,1:(end-1));
            xy0 = coords_full(1:2,end);
            xy = xy1 - xy0;
       end
       
       function BH = BH_from_xy_cal_0(ps, xy, BH_at_xyorigin)
            if nargin < 3
                xy0 = [0; 0];
            else
                xy0_full = ps.cal_0.meas.xy_out_T_BH_in *...
                      [BH_at_xyorigin; 1];
                xy0 = xy0_full(1:2,1);
            end
            Tmat = ps.cal_0.meas.bh_out_T_xy_in;
            coords_full = Tmat*[xy0+xy ; ones(1,size(xy,2))];
            BH = coords_full(1:2,:);
       end
       
       function [ BH_to_xy, xy_to_BH] = get_3pt_transform(ps)
           T = ps.(ps.current_set).meas;
           switch T.ang_units
                case 'deg'
                    ang_scale = pi/180;
                case 'rad'
                    ang_scale = 1;
                otherwise
                    error('Angular unit must be deg or rad');
           end
           
           % Get the transform from BH to xy
            % first set the origin point:
            A0 = [1 0 -T.BH_base_point(1);...
                  0 1 -T.BH_base_point(2);...
                  0 0 1];
            % then scale the BHx and BHy, and flip in y if required:
            flip_sign = (T.flip_y * -1) + (~T.flip_y * 1);
            A1 = [(T.BHx_measured_shift/T.BHx_current_shift) 0 0;...
                  0 (T.BHy_measured_shift/T.BHy_current_shift)*flip_sign 0;...
                  0 0 1];
            % Then shear the grid:
            shear_shift = tan( (pi/2) - (T.inc_theta*ang_scale) );
            A2 = [1 flip_sign*shear_shift 0;...
                  0 1 0
                  0 0 1];
            % Then scale in y, to match the shear angle:
            y_shear_scale = cos(shear_shift);
            A3 = [1 0             0;...
                  0 y_shear_scale 0;...
                  0 0             1];

            % Then rotate the grid:
            theta  = -T.x_theta*ang_scale;
            A4 = [cos(theta) sin(theta) 0;...
                 -sin(theta) cos(theta) 0;...
                 0 0 1];
            BH_to_xy = A4*A3*A2*A1*A0;
            xy_to_BH = inv(BH_to_xy);
           
       end
                     
       function find_offsets(ps, method_str, base_ind)
            if nargin < 3
                base_ind = 1;
            end
            if nargin < 2
                method_str = 'xcorr';
            end
            switch method_str
               case 'xcorr'
                   ps.(ps.current_set).base_ftRef = fft2(ps.(ps.current_set).ims{base_ind});
                   comp_method = @(im_in) ps.xcrr2(im_in);
               case 'sub-pix'
                   ps.(ps.current_set).base_ftRef = fft2(ps.(ps.current_set).ims{base_ind});                   
                   comp_method = @(im_in) ps.dftreg(im_in, 10);
               case 'sub-pix thresh'
                   % needs image processing toolbox
                  max_im = max(ps.(ps.current_set).ims{base_ind}(:)); 
                  thresh = graythresh(ps.(ps.current_set).ims{base_ind}/...
                                  max_im );                  
                  ps.(ps.current_set).base_ftRef = ...
                      fft2(double(ps.(ps.current_set).ims{base_ind} >...
                                     (thresh* max_im) ));
                  pass_wrapper = @(im_in) double( im_in> (thresh*max_im) );
                  comp_method = @(im_in) ps.dftreg(pass_wrapper(im_in), 10);
               otherwise
                   error('Unknown registration method');
           end

           ps.(ps.current_set).raw_offsets = ps.reg_shift(comp_method);             

       end
       
       function scale_offsets(ps)
            [scale_x, scale_y] = ps.(ps.current_set).tags{1}.TheData.Calibrations.Dimensions.Scale;
            ps.(ps.current_set).scaled_offsets.x = scale_x*ps.(ps.current_set).sorted_offsets.x; 
            ps.(ps.current_set).scaled_offsets.y = scale_y*ps.(ps.current_set).sorted_offsets.y;           
       end
        
       function offsets = reg_shift(ps, comp_method)
            num_im = length(ps.(ps.current_set).ims);
            offsets.x = zeros(1,num_im);
            offsets.y = zeros(1,num_im);
            for ii = 1:num_im             
                xy_shift = comp_method(ps.(ps.current_set).ims{ii});
                offsets.x(ii) = xy_shift.x;
                offsets.y(ii) = xy_shift.y;
            end
       end
        
       function sort_offsets(ps)
            
            raw_offsets = ps.(ps.current_set).raw_offsets;
            arr_size = size(ps.(ps.current_set).ims{1});
           
            dim1 = arr_size(1);
            dim2 = arr_size(2);
            num_pts = numel(raw_offsets.x);
            full_pts.x = zeros(1,4*num_pts);
            full_pts.y = zeros(1,4*num_pts);
            for ii = 1:length(raw_offsets.x)
                x_ur = raw_offsets.x(ii);
                y_ur = raw_offsets.y(ii);

                x_pts = [x_ur + dim2, x_ur + dim2, x_ur, x_ur];
                y_pts = [y_ur + dim1, y_ur, y_ur + dim1, y_ur];
                full_pts.x( (((ii-1)*4)+1):(ii*4)) = x_pts;
                full_pts.y( (((ii-1)*4)+1):(ii*4)) = y_pts;
            end
            f1_sel = figure; plot(full_pts.x, full_pts.y, 'x'); axis equal;
            disp('Select Contiguous Region');
            rect_region = getrect(f1_sel);
            close(f1_sel);
            % set signs for the quadrants:
            x_bnd = rect_region(1);
            y_bnd = rect_region(2);
            offsets.x = raw_offsets.x;
            offsets.y = raw_offsets.y;
            offsets.x(offsets.x >= x_bnd) = -(dim2 - offsets.x(offsets.x >= x_bnd));
            offsets.y(offsets.y >= y_bnd) = -(dim1 - offsets.y(offsets.y >= y_bnd));
            offsets.x = offsets.x - offsets.x(1);
            offsets.y = offsets.y - offsets.y(1);
            
            ps.(ps.current_set).sorted_offsets = offsets;

       end
       
       function offsets = xcrr2(ps, A)
            size_in = size(A,1);
            C = ifftshift(ifft2(fft2(A).*conj(ps.(ps.current_set).base_ftRef)));
            [~, ind] = max(C(:));
            [y_c, x_c ] = ind2sub(size(C), ind);
            offsets.x = x_c-((size_in/2)+1);
            offsets.y = -(y_c-((size_in/2)+1));     
       end
       
       function offsets = dftreg(ps, A, sub_pix)
           if nargin < 3
               sub_pix = 10;
           end
           a = dftregistration(ps.(ps.current_set).base_ftRef ,fft2(A), sub_pix); 
           offsets.x = -a(4);
           offsets.y =  a(3); 
       end
       
       function clear_load_list(ps, set_type)
           ps.(set_type).filelist_to_load = [];
       end
       
       function add_to_loadlist(ps, set_type ,newfilename, position, seq_number)
           if ~isempty(newfilename)
           base_inc = 500;
           if ~isfield(ps.(set_type),'filelist_to_load') || isempty(ps.(set_type).filelist_to_load)
               ps.(set_type).filelist_to_load.filenames = cell(base_inc, 1);
               ps.(set_type).filelist_to_load.loaded = zeros(base_inc,1);
               ps.(set_type).filelist_to_load.positions = zeros(2, base_inc);
               ps.(set_type).filelist_to_load.seq_number = zeros(1, base_inc);
               ps.(set_type).filelist_to_load.num_files_in_list = 0;
%                ps.(set_type).filelist_to_load.last_added = 0;
           end
           if ps.(set_type).filelist_to_load.num_files_in_list == ...
                   ps.(set_type).filelist_to_load.num_files_in_list + base_inc;
               ps.(set_type).filelist_to_load.filenames = ...
                   [ps.(set_type).filelist_to_load.filenames;...
                    cell(base_inc, 1)]; % grow the array
               ps.(set_type).filelist_to_load.loaded = ...
                   [ps.(set_type).filelist_to_load.loaded;...
                    zeros(base_inc,1)];
               ps.(set_type).filelist_to_load.positions = ...
                   [ps.(set_type).filelist_to_load.positions, zeros(2, base_inc)];
               ps.(set_type).filelist_to_load.seq_number = ...
                   [ps.(set_type).filelist_to_load.seq_number, zeros(1, base_inc)]; 
           end
           current_num = ps.(set_type).filelist_to_load.num_files_in_list + 1;
           % check it is not already in the file list
           if sum(strcmpi(newfilename,ps.(set_type).filelist_to_load.filenames)) == 0
               ps.(set_type).filelist_to_load.filenames{current_num} = newfilename;
               ps.(set_type).filelist_to_load.positions(:,current_num) = position;
               ps.(set_type).filelist_to_load.seq_number(1,current_num) = seq_number;
               ps.(set_type).filelist_to_load.num_files_in_list = 1 + ...
                   ps.(set_type).filelist_to_load.num_files_in_list;
               disp(['Added new file on to load list: ',newfilename]);
           else
               disp('Requested to loaded a file that is already on the load list');
           end
           end
       end
       
       function drift_scan_posns(ps, drift_mag_nm, drift_angle_deg, operate_on)
           if nargin < 4
               operate_on = 'scan_posns';
           end
           max_ind = size(ps.dps.scan_posns, 2);
           % drift vector is the vector of beam movement over the sample.
           % i.e. if image on screen appears to be moving down and to the
           % right, then effectively the beam is moving upwards and to the
           % left as time goes by. In this example for up and left the
           % position angle would be in the second quadrant, ie 90<angle<180
           drift_cx = (drift_mag_nm/max_ind) * exp(1j*2*pi*(drift_angle_deg/360)); 
           % Divide by max_ind to get per-position drift factor
           ps.dps.drifted_posns = ps.dps.(operate_on);
           for n = 1:max_ind
               ps.dps.drifted_posns(1,n) = ps.dps.(operate_on)(1,n) + real(drift_cx)*(n-1);
               ps.dps.drifted_posns(2,n) = ps.dps.(operate_on)(2,n) + imag(drift_cx)*(n-1);
           end

       end
       
       function [drift_mag_nm, drift_ang_deg] = calc_drift_vector(ps, pre_exp_im, post_exp_im, first_dp, last_dp)
           [~, drift_ang_deg, drift_mag_nmpersec] = ps.calc_drift(pre_exp_im, post_exp_im);
           tem_mode = false;
           if nargin == 5
               if ischar(first_dp) && ischar(last_dp)
                   if tem_mode
                       % im1 = importDM(first_dp);
                       % im2 = importDM(last_dp);
                       % date1 = im1.AcquisitionStartTimeStamp;
                       % date2 = im2.AcquisitionStartTimeStamp;
                       % elapsed_time = etime(datevec(date2), datevec(date1));
                       elapsed_time = 10; % just for testing purposes!
                       warning('Drift Calc currently not functional: need to put new DM3 reader that can extract time stamps! \n')

                   else
                       elapsed_time = first_dp; % just put the time here for now...
                   end
                   drift_mag_nm = drift_mag_nmpersec*elapsed_time;
               end
           end
           if nargin == 3
               % need to write a func to extract time stamps from tags
               date1 = ps.dps.AcquisitionStartTimeStamp{1};
               date2 = ps.dps.AcquisitionStartTimeStamp{end};
               elapsed_time = etime(datevec(date2), datevec(date1));
               drift_mag_nm = drift_mag_nmpersec*elapsed_time;
           end                      
       end
       
       function rotate_dp_posns(ps, operate_on)
           if nargin == 1
               operate_on = 'drifted_posns';
           end
           rot_angle = (ps.recon_dat.p.rot_angle_deg /180)*pi;  %Rotation seems to be closer to 140 degrees

           R = [cos(rot_angle) -sin(rot_angle);...
                sin(rot_angle)  cos(rot_angle)];
            
           xy_posns = ps.dps.(operate_on);
           ps.dps.rotated_posns = R*xy_posns;
       end
       
              
       function scale_dp_posns(ps, sf, operate_on)
           if nargin == 2
               operate_on = 'drifted_posns';
           end

           ps.dps.scaled_posns = sf*ps.dps.(operate_on);
       end
       
       function create_RC_coords(ps, operate_on)
           if nargin == 1
               operate_on = 'rotated_posns';
           end
           xy_posns = ps.dps.(operate_on)*1e-9; % as by default positions 
           % are given in nm, recon_pixel_size is in m.
           % Note for going into recon space 'x' is row ind, 'y' is col ind
           ps.dps.recon_space_posns = [-xy_posns(2,:);...
                                        xy_posns(1,:)] /... 
                                        ps.recon_dat.recon_pixel_size;
           ps.offset_posns_and_create_object_canvas;
           ps.dps.dpPosR = ps.dps.recon_space_posns(1,:); % Row
           ps.dps.dpPosC = ps.dps.recon_space_posns(2,:); % Col           
       end
       
       function offset_posns_and_create_object_canvas(ps)
            RC_posns_pix = ps.dps.recon_space_posns;
            % put everything into pixels
            min_R = min(RC_posns_pix(1,:)); max_R = max(RC_posns_pix(1,:));
            min_C = min(RC_posns_pix(2,:)); max_C = max(RC_posns_pix(2,:));
            % shift all to be based around 0 and add offset to put each pattern away
            % from the edge
%             obj_offset = ceil(1.5*(ps.recon_dat.p.dpSize/2));
            obj_offset = ceil(1.5*(ps.recon_dat.p.dpSize));
            
            ps.recon_dat.row_offset = obj_offset - min_R;
            ps.recon_dat.col_offset = obj_offset - min_C;
            RC_posns_pix(1,:) = RC_posns_pix(1,:) - min_R + obj_offset ;
            RC_posns_pix(2,:) = RC_posns_pix(2,:) - min_C + obj_offset ;

            range_R = ceil(max_R - min_R);
            range_C = ceil(max_C - min_C);
            ps.recon_dat.object_0 = ones(ceil([range_R  + 2*obj_offset,...
                                               range_C +  2*obj_offset]));
%             ps.recon_dat.object_0 = ps.recon_dat.object_0.'; %2020 04 20
            ps.dps.recon_space_posns = RC_posns_pix;
            ps.dps.posn_params.min_R = min_R;
            ps.dps.posn_params.max_R = max_R;
            ps.dps.posn_params.min_C = min_C;
            ps.dps.posn_params.max_R = max_R;
            ps.dps.posn_params.obj_offset = obj_offset;
       end
       
       function src_struct = setup_recon_params(ps, varargin)

            % look for default values:
            input_info = {...
                'probe_intensity', [],        'Image to make guess of probe intensity from',           '%s',      'INT_FILE';
                'BeamVoltage',200e3,          'Beam Energy (eV):',                                     '%-10.2f',  'V_BEAM';...
                'dpSize',1024,                'Diffraction Pattern Size (pixels):',                    '%-10.3f',  'DPSIZE';...
                'illum_diam',30e-9,           'Illumination Diam (m):',                                '%-10.2f',  'ILL_DIAM'; ...
                'angle_per_pix',0.06270*1e-3, 'Angle Per Pix at Det (rads):',                          '%-10.1f',  'SNR'; ...
                'central_disc_diam', 50,      'Approximate disc diameter in pixels in the DP):',      '%-16.5e',  'CS'; ...
                'rot_angle_deg', 0,           'Rotation Angle between DP and image',                   '%-16.5e',  'CS'; ...
                'calc_probe', 0,              'Rotation Angle between DP and image',                   '%-16.5e',  'CS'; ...
                'def_sign', -1,                'The defocus sign (conv/diverge illum)',                '%-16.5e',  'CC';
                'physical_cam_pix', 28e-6,    'Pyhysical Size of DP Pixel at Camera Plane',                   '%-16.5e',  'CS'}; 

            % angle_per_pix = 0.06270*1e-3; % from recent camera length calibration for 0.5 m length at 200 kv
            % illum_diam = 14e-9; % typical
            
            inputfields = input_info(:,1:2);
            param_names = [input_info(:,1), input_info(:,5:end)];

            read_dat = read_prog_input(inputfields, param_names, varargin{:});
            ps.recon_dat.p  = ps.parse_and_check_ip(read_dat.supplied_p, input_info);

            beam = RelatBeam(ps.recon_dat.p.BeamVoltage, 0);
            ps.recon_dat.lambda = ...
                beam.BeamLambdaRel;
            ps.recon_dat.beam_half_angle = ...
                ps.recon_dat.p.angle_per_pix * ps.recon_dat.p.central_disc_diam/2;
            ps.recon_dat.angular_range  = ...
                ps.recon_dat.p.dpSize*ps.recon_dat.p.angle_per_pix;
            ps.recon_dat.recon_pixel_size = ...
                ps.recon_dat.lambda/(ps.recon_dat.angular_range );

            % approximate defocus:
            ps.recon_dat.defocus = ps.recon_dat.p.def_sign * ...
                (ps.recon_dat.p.illum_diam/2)/ps.recon_dat.beam_half_angle;
            
            if ps.recon_dat.p.calc_probe
                [probe_pf] = ps.perfect_probe(...
                                   ps.recon_dat.p.dpSize, ...
                                   ps.recon_dat.recon_pixel_size, ...
                                   ps.recon_dat.defocus, ...
                                   ps.recon_dat.beam_half_angle, ...
                                   ps.recon_dat.lambda);

                
                if isempty(ps.recon_dat.p.probe_intensity)
                    ps.recon_dat.probe_pf = probe_pf;
                else
                    probe_modulus = ps.probe_guess_from_dm3(...
                        ps.recon_dat.p.probe_intensity, 'dm3', ...
                        ps.recon_dat.recon_pixel_size, ps.recon_dat.p.dpSize, ps.recon_dat.p.rot_angle_deg);
                    probe_angle = angle(probe_pf);
                    ps.recon_dat.probe_pf = probe_modulus .* exp(1i*probe_angle);
                end
                
                ps.recon_dat.probe_pf(isnan(ps.recon_dat.probe_pf)) = 0; 
            end 
            src_struct = ps.recon_dat.p;
       end
       
      
       
       function save_probe_pf(ps, savelocation)
           ps.recon_dat.ptychoshelves.initial_probe_file = savelocation;
           p.binning = 0;
           p.detector.binning = 0;           
           probe = ps.recon_dat.probe_pf;
           save(savelocation, 'p','probe');
       end
       
       
       function [p, eng] = conv_ptychoset_cSAXS_MSlc(ps, set_params, set_number, recipe, recipe_num, local_exec_params)
           
            if nargin < 6
                local_exec_params.verbose_level = 2; 
            end
           
            if nargin < 5
                recipe_num = 1;
            end
            
            simple_check_recipe_field = @(r_num, fn) ...
                isfield(recipe.params(r_num), fn) && ...
                ~isnan(recipe.params(r_num).(fn)) && ...
                (recipe.params(r_num).(fn)) > 0;
            
            Ndpx = ps.recon_dat.p.dpSize;
            HT =  ps.recon_dat.p.BeamVoltage/1e3;
            d_alpha =  ps.recon_dat.p.angle_per_pix;
            idx_scan=1; % scan area index
            
            p = struct();
            p.   verbose_level = local_exec_params.verbose_level;          % verbosity for standard output (0-1 for loops, 2-3 for testing and adjustments, >= 4 for debugging)
            p.   use_display = false;        % global switch for display, if [] then true for verbose > 1
            p.   scan_number = [idx_scan];   % Multiple scan numbers for shared scans
            
            % Geometry
            p.   z = 1/d_alpha ;             % 1/delta_angle in rad^-1, for electron ptychography 
            p.   asize = [Ndpx,Ndpx];        % Diffr. patt. array size
            p.   ctr = [fix(Ndpx/2)+1, fix(Ndpx/2)+1];    % Diffr. patt. center coordinates (y,x) (empty means middle of the array); e.g. [100 207;100+20 207+10];
            p.   prop_regime = 'farfield';                % propagation regime: nearfield, farfield (default), !! nearfield is supported only by GPU engines 
            p.   energy = HT ;                            % Energy (in keV), leave empty to use spec entry mokev
            p.   electron = true;                         % for electron ptychography, added by ZC
            
            p.   affine_matrix = [1 , 0; 0, 1] ; % Applies affine transformation (e.g. rotation, stretching) to the positions (ignore by = []). Convention [yn;xn] = M*[y;x]. For flOMNI we found in June 2019: = [1 , 0.0003583 ; 5.811e-05 , 1 ]; for OMNY we found in October 2018: = [1 0;tan(0.4*pi/180) 1]; laMNI in June 2018  [1,0.0154;-0.0017,1.01]; laMNI in August [1.01 0.0031; -0.0018 1.00] 
            % Scan meta data
            p.   src_metadata = 'none';          % source of the meta data, following options are supported: 'spec', 'none' , 'artificial' - or add new to +scan/+meta/
            p.   queue.lockfile = false;         % If true writes a lock file, if lock file exists skips recontruction
            
            % Data preparation
            p.   detector.name = 'empad';                 % see +detectors/ folder 
            p.   detector.check_2_detpos = [];            % = []; (ignores)   = 270; compares to dettrx to see if p.ctr should be reversed (for OMNY shared scans 1221122), make equal to the middle point of dettrx between the 2 detector positions
            p.   detector.data_prefix = '';               % Default using current eaccount e.g. e14169_1_
            p.   detector.binning = false;                % = true to perform 2x2 binning of detector pixels, for binning = N do 2^Nx2^N binning
            p.   detector.upsampling = false;             % upsample the measured data by 2^data_upsampling, (transposed operator to the binning), it can be used for superresolution in nearfield ptychography or to account for undersampling in a far-field dataset
            p.   detector.burst_frames = 1;               % number of frames collected per scan position
            p.   detector.dp_dat_filename = set_params.file_name{set_number};
            
            p.   prepare.data_preparator = 'matlab_aps';  % data preparator; 'python' or 'matlab' or 'matlab_aps'
            p.   prepare.auto_prepare_data = true;        % if true: prepare dataset from raw measurements if the prepared data does not exist
            p.   prepare.force_preparation_data = true;   % Prepare dataset even if it exists, it will overwrite the file % Default: @prepare_data_2d
            p.   prepare.store_prepared_data = false;     % store the loaded data to h5 even for non-external engines (i.e. other than c_solver)
            p.   prepare.prepare_data_function = '';      % (used only if data should be prepared) custom data preparation function handle;
            p.   prepare.auto_center_data = false;        % if matlab data preparator is used, try to automatically center the diffraction pattern to keep center of mass in center of diffraction
            
            % Scan positions
            p.   src_positions = 'hdf5_pos';          % 'spec', 'orchestra', 'load_from_file', 'matlab_pos' (scan params are defined below) or add new position loaders to +scan/+positions/
            p.   positions_file = '';    %Filename pattern for position files, Example: ['../../specES1/scan_positions/scan_%05d.dat']; (the scan number will be automatically filled in)

            
            % scan parameters for option src_positions = 'matlab_pos';
            scan_string_format = '%02d';
            p.   scan.type = 'custom';                     % {'round', 'raster', 'round_roi', 'custom'}
            p.   scan.roi_label = [];                      % For APS data
            p.   scan.format = scan_string_format;         % For APS data format for scan directory generation
            p.   scan.custom_positions_source = ps.recon_dat.ptychoshelves.crds_file;   % custom: a string name of a function that defines the positions; also accepts mat file with entry 'pos', see +scans/+positions/+mat_pos.m
            p.   scan.custom_params = [];                  % custom: the parameters to feed to the custom position function.
            p.   scan.trans_coords = set_params.trans_coords(set_number); % parameter to allow transformation of position coordinates.
            p.   scan.rotation = set_params.meas_rot_est(set_number) - set_params.rot_applied(set_number);
            p.   scan.scale_in_ps = set_params.scale_in_ps(set_number);
            
            % I/O
            p.   prefix = local_exec_params.recon_savename;    % For automatic output filenames. If empty: scan number;
            p.   suffix = '';                             % Optional suffix for reconstruction 
            p.   scan_string_format = scan_string_format; % format for scan string generation, it is used e.g for plotting and data saving 
            if any(strcmpi(set_params.Properties.VariableNames, 'recon_savepath')) && ~isempty(set_params.recon_savepath{set_number}) && ischar(set_params.recon_savepath{set_number})
                reconsavepath = set_params.recon_savepath{set_number};
            else
                reconsavepath = set_params.base_dir{set_number};
            end
                        
            p.   base_path = set_params.base_dir{set_number};   % base path : used for automatic generation of other paths
            
            p.   specfile = '';                                         % Name of spec file to get motor positions and check end of scan, defaut is p.spec_file == p.base_path;
            if isfield(local_exec_params,'ptycho_matlab_path')
                p.   ptycho_matlab_path = local_exec_params.ptycho_matlab_path;
            else
                p.   ptycho_matlab_path = 'U:\Ptycho\PtychoShelves_EM-1.0\ptycho';  % cSAXS ptycho package path
            end
            if isfield(local_exec_params,'cSAXS_matlab_path')
                p.   cSAXS_matlab_path = local_exec_params.cSAXS_matlab_path;
            else
                p.   cSAXS_matlab_path = 'U:\Ptycho\PtychoShelves_EM-1.0';  % cSAXS base package path    
            end
                                    
            p.   raw_data_path{1} = '';                                 % Default using compile_x12sa_filename, used only if data should be prepared automatically
            p.   prepare_data_path = fullfile(reconsavepath,'intermediates');    % Default: base_path + 'analysis'. Other example: '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/'; also supports %u to insert the scan number at a later point (e.g. '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/S%.5u')
            p.   prepare_data_filename = [];                            % Leave empty for default file name generation, otherwise use [sprintf('S%05d_data_%03dx%03d',p.scan_number(1), p.asize(1), p.asize(2)) p.prep_data_suffix '.h5'] as default 
            p.   save_path{1} = fullfile(reconsavepath,'intermediates');         % Default: base_path + 'analysis'. Other example: '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/'; also supports %u to insert the scan number at a later point (e.g. '/afs/psi.ch/project/CDI/cSAXS_project/analysis2/S%.5u')
            p.   io.default_mask_file = '';                             % load detector mask defined in this file instead of the mask in the detector packages, (used only if data should be prepared) 
            p.   io.default_mask_type = 'binary';                       % (used only if data should be prepared) ['binary', 'indices']. Default: 'binary' 
            p.   io.file_compression = 0;                               % reconstruction file compression for HDF5 files; 0 for no compression
            p.   io.data_compression = 3;                               % prepared data file compression for HDF5 files; 0 for no compression
            p.   io.load_prep_pos = false;                              % load positions from prepared data file and ignore positions provided by metadata
           
            % Reconstruction
            % Initial iterate object
            p.   model_object = ps.recon_dat.ptychoshelves.model_object;       % Use model object, if false load it from file % ***------<<<<<<<<
            p.   model.object_type = ps.recon_dat.ptychoshelves.object_type;   % specify how the object shall be created; use 'rand' for a random initial guess; use 'amplitude' for an initial guess based on the prepared data
            p.   initial_iterate_object_file{1} = ps.recon_dat.ptychoshelves.initial_iterate_object_file{1}; % ***------<<<<<<<<
                   % use this mat-file as initial guess of object, it is possible to use wild characters and pattern filling, example: '../analysis/S%05i/wrap_*_1024x1024_1_recons*'
                   % PtychoProcessFromTable sets this as 'from_ptychoset_struct' in most cases, or as '' when rand is used.
            p.   distrib_objphase = ps.recon_dat.ptychoshelves.distrib_objphase;
            p.   mean_filter_obj = ps.recon_dat.ptychoshelves.mean_filter_obj;
            p.   tie_unwrap_obj = ps.recon_dat.ptychoshelves.tie_unwrap_obj;
            p.   compute_cumu_phases = ps.recon_dat.ptychoshelves.compute_cumu_phases;
            p.   cumu_phase_regulation = ps.recon_dat.ptychoshelves.cumu_phase_regulation;
            if simple_check_recipe_field(recipe_num, 'interleave_obj_phase')
                p.   interleave_obj_phase = recipe.params(recipe_num).interleave_obj_phase;
            else
                p.   interleave_obj_phase = 0;
            end

            % Initial iterate probe
            p.   model_probe = false;                                   % Use model probe, if false load it from file
            % My assumption (to be checked by looking through code /
            % experimenting) is that the group below are not used as the
            % above is false. Thus the group below is just a place holder
            % for now.
            p.   model.probe_is_focused = true;                         % Model probe is focused (false: just a pinhole)
            p.   model.probe_central_stop = true;                       % Model central stop
            p.   model.probe_diameter = 170e-6;                         % Model probe pupil diameter
            p.   model.probe_central_stop_diameter = 50e-6;             % Model central stop diameter
            p.   model.probe_zone_plate_diameter = 170e-6;              % Model probe zone plate diameter
            p.   model.probe_outer_zone_width = [];                     % Model probe zone plate outermost zone width (not used if not a focused probe) 
            p.   model.probe_propagation_dist = 3e-3;                   % Model probe propagation distance (pinhole <-> sample for unfocused, focal-plane <-> sample for focused)
            p.   model.probe_focal_length = 51e-3;                      % Model probe focal length (used only if model_is_focused is true
                                                                        %   AND model_outer_zone_width is empty)
            p.   model.probe_upsample = 10;                             % Model probe upsample factor (for focused probes)
            
            %Use probe from this mat-file (not used if model_probe is true)
%             p.   initial_probe_file = fullfile(p.base_path,sprintf(p.scan.format, p.scan_number),'probe_initial.mat');
            p.   initial_probe_file = ps.recon_dat.ptychoshelves.initial_probe_file; % ***---------------------------------------------------------------------<<<<<<<<<
            p.   probe_file_propagation = 0.0e-3;                            % Distance for propagating the probe from file in meters, = 0 to ignore
            
            % Shared scans - Currently working only for sharing probe and object
            p.   share_probe  = 0;                                      % Share probe between scans. Can be either a number/boolean or a list of numbers, specifying the probe index; e.g. [1 2 2] to share the probes between the second and third scan. 
            p.   share_object = 0;                                      % Share object between scans. Can be either a number/boolean or a list of numbers, specifying the object index; e.g. [1 2 2] to share the objects between the second and third scan. 

            % Modes
            p.   probe_modes  = recipe.params(recipe_num).N_probes;       % Number of coherent modes for probe
            p.   object_modes = 1;                                      % Number of coherent modes for object
            % Set pixels in the dp below this level to zero, i.e. assume
            % that below this level it is incoherent noise.
            p.   dp_Pedestals_incoherent = recipe.params(recipe_num).dp_Pedestals_incoherent;           
            
            if simple_check_recipe_field(recipe_num, 'edgezeros')
                p.   edgezeros = recipe.params(recipe_num).edgezeros;
            else
                p.   edgezeros = 0;
            end
            
            if any(strcmpi(fieldnames(recipe.params), 'border_mask'))
                if isa(recipe.params(recipe_num).border_mask,'string')
                    if strcmpi(recipe.params(recipe_num).border_mask, 'imdef')
                        p.   border_mask = 'imdef'; % this sets it as char type. Before this it would be string type.
                    else
                        p.   border_mask = str2num(recipe.params(recipe_num).border_mask);
                    end
                else
                    p.   border_mask = recipe.params(recipe_num).border_mask;
                end
            else
                p.   border_mask = 0;
            end
            
            if any(strcmpi(fieldnames(recipe.params), 'cross_mask'))
                if isa(recipe.params(recipe_num).cross_mask,'string')
                    p.   cross_mask = str2num(recipe.params(recipe_num).cross_mask);
                else
                    p.   cross_mask = recipe.params(recipe_num).cross_mask;
                end
            else
                p.   cross_mask = 0;
            end
            
            % Mode starting guess
            p.   mode_start_pow = 0.02;                               % Normalized intensity on probe modes > 1. Can be a number (all higher modes equal) or a vector
            p.   mode_start = 'herm';                                   % (for probe) = 'rand', = 'herm' (Hermitian-like base), = 'hermver' (vertical modes only), = 'hermhor' (horizontal modes only)
            p.   ortho_probes = true;                                   % orthogonalize probes after each engine
            p.   object_regular = 0;                                    % should be smaller than 1/8, smooth object amplitude, usefull for many layers
                                                                        % NOTE: object_regular = reg_mu (in load_from_p.m line 102) -> 
                                                                        % applies smoothness constraint on the object amplitude, apply_smoothness_constraint(object, object_regular)
                                                                        % i.e. reg_mu (set below also controls object_regular)
            
            %% Plot, save and analyze
            p.   plot.prepared_data = false;                         % plot prepared data
            p.   use_display = false;
            p.   plot.windowautopos = true;
            p.   plot.show_layers = 1;
            p.   plot.show_layers_stack = true;
            p.   plot.conjugate = false;
            p.   plot.horz_fact = 2;
            p.   plot.realaxes = true;
            p.   plot.fov_box = true;
            p.   plot.fov_box_color = 'r';
            p.   plot.positions = true;
            p.   save.external = true;                             % Use a new Matlab session to run save final figures (saves ~6s per reconstruction). Please be aware that this might lead to an accumulation of Matlab sessions if your single reconstruction is very fast.
            p.   save.store_images = false;                              % Write preview images containing the final reconstructions in [p.base_path,'analysis/online/ptycho/'] if p.use_display = 0 then the figures are opened invisible in order to create the nice layout. It writes images in analysis/online/ptycho
            p.   save.store_images_intermediate = false;                % save images to disk after each engine
            p.   save.store_images_ids = 1:6;                           % identifiers  of the figure to be stored, 1=obj. amplitude, 2=obj. phase, 3=probes, 4=errors, 5=probes spectrum, 6=object spectrum
            p.   save.store_images_format = 'png';                      % data type of the stored images jpg or png 
            p.   save.store_images_dpi = 150;                           % DPI of the stored bitmap images 
%             p.   save.exclude = {'fmag', 'fmask', 'illum_sum'};         % exclude variables to reduce the file size on disk
            p.   save.exclude = {'fmag', 'fmask'};         % exclude variables to reduce the file size on disk
            p.   save.save_reconstructions_intermediate = false;        % save final object and probes after each engine
            p.   save.save_reconstructions = false;                      % save reconstructed object and probe when full reconstruction is finished 
            p.   save.output_file = 'h5';                               % data type of reconstruction file; 'h5' or 'mat'

% ---- The Engine parameters ----
                % --------- GPU engines  -------------   See for more details: OdstrÄ?il M, et al., Optics express. 2018 Feb 5;26(3):3108-23.                                                  
            eng = struct();                        % reset settings for this engine 
            eng. name = 'GPU';       
            eng. use_gpu = true;                   % if false, run CPU code, but it will get very slow 
            
            eng. keep_on_gpu = true;              % keep data + projections on GPU, false is useful for large data if DM is used
            eng. compress_data = false;            % use automatic online memory compression to limit need of GPU memory
            eng. gpu_id = [];                      % default GPU id, [] means choosen by matlab
            eng. check_gpu_load = true;            % check available GPU memory before starting GPU engines 

            % general 
            eng. number_iterations = recipe.params(recipe_num).iterations;          % number of iterations for selected method
            if isfinite(recipe.params(recipe_num).Np_presolve)
                eng. asize_presolve = recipe.params(recipe_num).Np_presolve;
            else
                eng. asize_presolve = ps.recon_dat.p.dpSize;                            % CHECK what 'asize_presolve' actaully does
            end

            eng. method = recipe.params(recipe_num).algorithm;                       % choose GPU solver: DM, ePIE, hPIE, MLc, Mls, -- recommended are MLc and MLs
            eng. opt_errmetric = 'L1';            % optimization likelihood - poisson, L1
            eng. grouping = recipe.params(recipe_num).para_grouping;        % size of processed blocks, larger blocks need more memory but they use GPU more effeciently, !!! grouping == inf means use as large as possible to fit into memory 
                                                   % * for hPIE, ePIE, MLs methods smaller blocks lead to faster convergence, 
                                                   % * for MLc the convergence is similar 
                                                   % * for DM is has no effect on convergence
            eng. probe_modes  = p.probe_modes;                % Number of coherent modes for probe
            eng. object_change_start = recipe.params(recipe_num).object_reconstruct_start;          % Start updating object at this iteration number
            eng. probe_change_start = recipe.params(recipe_num).probe_reconstruct_start;          % Start updating probe at this iteration number
            

            % regularizations
            if simple_check_recipe_field(recipe_num, 'reg_mu')
                eng. reg_mu = recipe.params(recipe_num).reg_mu;
            else
                eng. reg_mu = 0;                       % Regularization (smooting) constant ( reg_mu = 0 for no regularization)
            end
            if simple_check_recipe_field(recipe_num, 'delta')
                eng. delta = recipe.params(recipe_num).delta;
            else
                eng. delta = 0;                        % press values to zero out of the illumination area in th object, usually 1e-2 is enough
            end
            eng. positivity_constraint_object = 0; % enforce weak (relaxed) positivity in object, ie O = O*(1-a)+a*|O|, usually a=1e-2 is already enough. Useful in conbination with OPRP or probe_fourier_shift_search  

            eng. apply_multimodal_update = false; % apply all incoherent modes to object, it can cause isses if the modes collect some crap 
            eng. probe_backpropagate = 0;         % backpropagation distance the probe mask, 0 == apply in the object plane. Useful for pinhole imaging where the support can be applied  at the pinhole plane
            if simple_check_recipe_field(recipe_num, 'probe_support_radius') 
                eng. probe_support_radius = recipe.params(recipe_num).probe_support_radius;       % Normalized radius of circular support, = 1 for radius touching the window    
            else
                eng. probe_support_radius = [];
            end
            eng. probe_support_fft = false;       % assume that there is not illumination intensity out of the central FZP cone and enforce this contraint. Useful for imaging with focusing optics. Helps to remove issues from the gaps between detector modules.
            eng. probe_support_tem = false;        % assume a binary mask aperture for TEM, generated from initial probe, by Zhen Chen
            if simple_check_recipe_field(recipe_num, 'probe_support_apt_radius')
                eng. probe_support_apt_radius = recipe.params(recipe_num).probe_support_apt_radius;       % assume a binary mask aperture for TEM, generated from initial probe, by Zhen Chen
            else
                recipe.params(recipe_num).probe_support_apt_radius = Inf;
            end
            % basic recontruction parameters 
            % PIE / ML methods                    % See for more details: OdstrÄ?il M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
            eng. beta_object = recipe.params(recipe_num).rp_Object_updateFactor;                 % object step size, larger == faster convergence, smaller == more robust, should not exceed 1
            eng. beta_probe = recipe.params(recipe_num).rp_Probe_updateFactor;                  % probe step size, larger == faster convergence, smaller == more robust, should not exceed 1
            eng. delta_p = recipe.params(recipe_num).lsq_delta_p;                   % LSQ dumping constant, 0 == no preconditioner, 0.1 is usually safe, Preconditioner accelerates convergence and ML methods become approximations of the second order solvers 
            eng. momentum = recipe.params(recipe_num).momentum;                    % add momentum acceleration term to the MLc method, useful if the probe guess is very poor or for acceleration of multilayer solver, but it is quite computationally expensive to be used in conventional ptycho without any refinement. 
            eng. beta_LSQ = recipe.params(recipe_num).beta_LSQ;                                      % The momentum method works usually well even with the accelerated_gradients option.  eng.momentum = multiplication gain for velocity, eng.momentum == 0 -> no acceleration, eng.momentum == 0.5 is a good value
                                                  % momentum is enabled only when par.Niter < par.accelerated_gradients_start;
            eng. accelerated_gradients_start = inf; % iteration number from which the Nesterov gradient acceleration should be applied, this option is supported only for MLc method. It is very computationally cheap way of convergence acceleration. 

            % ADVANCED OPTIONS                     See for more details: OdstrÄ?il M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
            % position refinement 
%             eng. apply_subpix_shift = true;       % apply FFT-based subpixel shift, it is automatically allowed for position refinement
            
            eng. apply_subpix_shift = recipe.params(recipe_num).rp_subPixel;   % apply FFT-based subpixel shift, important for good position refinement but it is slow
            if ~simple_check_recipe_field(recipe_num, 'rp_Pos_Random_enabled')
%             if ~recipe.params(recipe_num).rp_Pos_Random_enabled
               eng.probe_position_search = inf;       % iteration number from which position correction is started
            else
               if isnumeric(recipe.params(recipe_num).rp_Pos_Random_enabled)
                   eng.probe_position_search = recipe.params(recipe_num).rp_Pos_Random_enabled;
                   fprintf('Applying Position correction after iteration %.0f\n',recipe.params(recipe_num).rp_Pos_Random_enabled);
               end
            end

            eng. probe_geometry_model = {};  % list of free parameters in the geometry model, choose from: {'scale', 'asymmetry', 'rotation', 'shear'}
            eng. probe_position_error_max = inf; % in meters, maximal expected random position errors, probe prositions are confined in a circle with radius defined by probe_position_error_max and with center defined by original positions scaled by probe_geometry_model
            eng. apply_relaxed_position_constraint = false;  % added by YJ
            % multilayer extension 
            if recipe.params(recipe_num).use_multislice
                delta_z = recipe.params(recipe_num).delta_z_angstrom*1e-10;
                Nlayers = recipe.params(recipe_num).Nlayers;
                eng. delta_z = delta_z*ones(Nlayers,1);             % if not empty, use multilayer ptycho extension , see ML_MS code for example of use, [] == common single layer ptychography , note that delta_z provides only relative propagation distance from the previous layer, ie delta_z can be either positive or negative. If preshift_ML_probe == false, the first layer is defined by position of initial probe plane. It is useful to use eng.momentum for convergence acceleration 
                eng. regularize_layers = recipe.params(recipe_num).regularize_layers;            % multilayer extension: 0<R<<1 -> apply regularization on the reconstructed object layers, 0 == no regularization, 0.01 == weak regularization that will slowly symmetrize information content between layers 
%                 eng. preshift_ML_probe = false;         % multilayer extension: if true, assume that the provided probe is reconstructed in center of the sample and the layers are centered around this position             
                eng. preshift_ML_probe = true;         % multilayer extension: if true, assume that the pr
            else
                eng.delta_z = [];
            end
            % other extensions - backgrounding models
            %%
            if simple_check_recipe_field(recipe_num, 'background_detection')
                eng. background_detection = recipe.params(recipe_num).background_detection; % the iteration number at which packground_detection starts to be applied. To turn off set as Inf.
            else
                eng. background_detection = inf;
            end
            if simple_check_recipe_field(recipe_num, 'background')
                eng. background = recipe.params(recipe_num).background;
            else
                eng. background = 0;               % average background scattering level, for OMNI values around 0.3 for 100ms, for flOMNI <0.1 per 100ms exposure, see for more details: Odstrcil, M., et al., Optics letters 40.23 (2015): 5574-5577.
            end
            %%
            if simple_check_recipe_field(recipe_num, 'background_width')
                eng. background_width = recipe.params(recipe_num).background_width;            % width of the background function in pixels,  inf == flat background, background function is then convolved with the average diffraction pattern in order to account for beam diversion 
            else
                eng. background_width = inf;
            end
            
            if simple_check_recipe_field(recipe_num, 'diff_pattern_blur')
                eng. diff_pattern_blur = recipe.params(recipe_num).diff_pattern_blur;            % width of the background function in pixels,  inf == flat background, background function is then convolved with the average diffraction pattern in order to account for beam diversion 
            else
                eng. diff_pattern_blur = 0;
            end
            
            
            
            
            eng. clean_residua = false;            % remove phase residua from reconstruction by iterative unwrapping, it will result in low spatial freq. artefacts -> object can be used as an residua-free initial guess for netx engine

            
%             eng. background = recipe.params(recipe_num).background;               % average background scattering level, for OMNI values around 0.3 for 100ms, for flOMNI <0.1 per 100ms exposure, see for more details: Odstrcil, M., et al., Optics letters 40.23 (2015): 5574-5577.
%             eng. background_width = recipe.params(recipe_num).background_width;           % width of the background function in pixels,  inf == flat background, background function is then convolved with the average diffraction pattern in order to account for beam diversion 
%             eng. clean_residua = false;            % remove phase residua from reconstruction by iterative unwrapping, it will result in low spatial freq. artefacts -> object can be used as an residua-free initial guess for netx engine

            % wavefront & camera geometry refinement     See for more details: OdstrÄ?il M, et al., Optics express. 2018 Feb 5;26(3):3108-23.
            if simple_check_recipe_field(recipe_num, 'probe_fourier_shift_search') 
                eng. probe_fourier_shift_search = recipe.params(recipe_num).probe_fourier_shift_search; % iteration number from which the engine will: refine farfield position of the beam (ie angle) from iteration == probe_fourier_shift_search
            else
                eng. probe_fourier_shift_search = inf;
            end
            eng. estimate_NF_distance = inf;       % iteration number from which the engine will: try to estimate the nearfield propagation distance using gradient descent optimization  
            if simple_check_recipe_field(recipe_num, 'detector_rotation_search') 
                eng. detector_rotation_search = recipe.params(recipe_num).detector_rotation_search;   % iteration number from which the engine will: search for optimal detector rotation, preferably use with option mirror_scan = true , rotation of the detector axis with respect to the sample axis, similar as rotation option in the position refinement geometry model but works also for 0/180deg rotation shared scans 
            else
                eng. detector_rotation_search = inf;
            end
            if simple_check_recipe_field(recipe_num,'detector_scale_search')  
                eng. detector_scale_search = recipe.params(recipe_num).detector_scale_search;      % iteration number from which the engine will: refine pixel scale of the detector, can be used to refine propagation distance in ptycho 
            else
                eng. detector_scale_search = inf; 
            end
            
            eng. variable_probe = false;           % Use SVD to account for variable illumination during a single (coupled) scan, see for more details:  Odstrcil, M. et al. Optics express 24.8 (2016): 8360-8369.
            eng. variable_probe_modes = 1;         % OPRP settings , number of SVD modes using to describe the probe evolution. 
            eng. variable_probe_smooth = 0;        % OPRP settings , enforce of smooth evolution of the OPRP modes -> N is order of polynomial fit used for smoothing, 0 == do not apply any smoothing. Smoothing is useful if only a smooth drift is assumed during the ptycho acquisition 
            eng. variable_intensity = false;       % account to changes in probe intensity

            % extra analysis
            eng. get_fsc_score = false;            % measure evolution of the Fourier ring correlation during convergence 
            eng. mirror_objects = false;           % mirror objects, useful for 0/180deg scan sharing -> geometry refinement for tomography, works only if 2 scans are provided 

            % custom data adjustments, useful for offaxis ptychography
            eng.auto_center_data = false;           % autoestimate the center of mass from data and shift the diffraction patterns so that the average center of mass corresponds to center of mass of the provided probe 
            eng.auto_center_probe = false;          % center the probe position in real space before reconstruction is started 
%             eng.custom_data_flip = [0,0,0];         % apply custom flip of the data [fliplr, flipud, transpose]  - can be used for quick testing of reconstruction with various flips or for reflection ptychography
            eng.custom_data_flip = [recipe.params(recipe_num).dp_FlipHorizontal, recipe.params(recipe_num).dp_FlipVertical, recipe.params(recipe_num).dp_Transpose];
            eng.apply_tilted_plane_correction = ''; % if any(p.sample_rotation_angles([1,2]) ~= 0),  this option will apply tilted plane correction. (a) 'diffraction' apply correction into the data, note that it is valid only for "low NA" illumination  Gardner, D. et al., Optics express 20.17 (2012): 19050-19059. (b) 'propagation' - use tilted plane propagation, (c) '' - will not apply any correction 

            eng.plot_results_every = inf;
            eng.save_results_every = recipe.params(recipe_num).save_results_every;
            eng.save_phase_image = true;
            eng.save_probe_mag = true;

%             resultDir = strcat(p.base_path,sprintf(p.scan.format, p.scan_number));
            resultDir = p.save_path{1};
            strcustom = local_exec_params.recon_savename;
            eng.fout =  fullfile(resultDir, strcustom);
            if ~exist(resultDir,'dir')
                fprintf('Making results directory: %s\n',resultDir);     
                mkdir(resultDir);
            end
                    
%{            
            % from previous version of LSQML
            if isa(ps.dps.ims_proc, 'uint16')
                inputs.diffraction = single(ps.dps.ims_proc);
            else
                inputs.diffraction = ps.dps.ims_proc;
            end
            if isfield(ps.dps,'masks')
                inputs.masks = ~ps.dps.masks;
            end
            if isfield(ps.dps,'edgezeros')
                inputs.edgezeros =  ps.dps.edgezeros;
            end
            
            inputs.object_orig = ps.recon_dat.object_0.'; % testing to see if this resolves the subregion problem with using Matlab ePIE implementation.
            inputs.probe_orig = ps.recon_dat.probe_pf.'; % testing to see if this resolves the subregion problem with using Matlab ePIE implementation.
            inputs.probe_positions_0 = [ps.dps.dpPosR',ps.dps.dpPosC'];
            inputs.probe_positions_0 = inputs.probe_positions_0  - mean(inputs.probe_positions_0);
            inputs.reconstruct_ind{1} = 1:length(inputs.probe_positions_0);
            inputs.Np_p = size(inputs.probe_orig);
            inputs.Np_o = size(inputs.object_orig);
            inputs.Npos = length(inputs.probe_positions_0);
            inputs.object{1} = inputs.object_orig;
            inputs.probe{1} = inputs.probe_orig;
            inputs.probe_positions = [];
            inputs.probe_support = [];
            inputs.modes = [];
            inputs.probe_evolution = [];            
%}
       end
       
    % end of methods   
    end
    
    methods (Static)
        
        function DP = maskDP(DP, DPmask)
            DP(~DPmask) = 0;
        end
        
        function param = conv_ptychosetrecipe_lsmql( recipe, iter)
           param.beta_object = recipe.params(iter).rp_Object_updateFactor;
           param.beta_probe = recipe.params(iter).rp_Probe_updateFactor;
           if recipe.params(iter).rp_Probe_update
               param.beta_probe = recipe.params(iter).rp_Probe_updateFactor;
           else
               param.beta_probe = 0;
           end
           param.Niter = recipe.iterations(iter);
           % 'defaults':
           param.probe_support_radius = [];  % radius of circular probe suppport, [] = none, useful for DM code 
           param.probe_support_distance = 0; % distance of the probe support from the object plane, 0 = none, inf = farfield
           param.Nprobes = 1;
           param.method = 'ePIE';
           param.grouping = 1;
           param.probe_reconstruct = 1;    % iteration when the probe reconstruction is started
           param.object_reconstruct = recipe.params(iter).object_reconstruct_start;   % iteration when the object reconstruction is started
           param.verbose_level = 1;        % verbosity of the solver , 0 = quiet, 4 = all
           %% DM
           param.pfft_relaxation = 0.05;  % relaxed modulus constraint for DM solver 
           param.probe_inertia = 0.1;
           %% ADVANCED OPTIONS       
           param.variable_probe = false;      % Use SVD to account for variable illumination during a single scan (Orthogonal probe relaxation)
           param.variable_probe_modes = 1;    % Number of SVD modes used for the orthogonal probe relaxation, ML method allows only one variable mode
           param.variable_intensity = false;  % Account for variable illumination intensity during a single scan

           param.apply_subpix_shift = recipe.params(iter).rp_subPixel;   % apply FFT-based subpixel shift, important for good position refinement but it is slow
           if ~recipe.params(iter).rp_Pos_Random_enabled
               param.probe_pos_search = inf;       % iteration number from which position correction is started
           else
               if isnumeric(recipe.params(iter).rp_Pos_Random_enabled)
                   param.probe_pos_search = recipe.params(iter).rp_Pos_Random_enabled;
                   fprintf('Applying Position correction after iteration %.0f',recipe.params(iter).rp_Pos_Random_enabled);
               end
           end

           param.apply_multimodal_update = false; % apply higher Thibault's modes to the object update, may result in lower object quality for real datasets
           param.beta_LSQ = true;   % use LSQ preditive step for ML method 
           param.delta_p = 0.1;     % LSQ damping constant, used only for MLs method                      
        end
       
        
        function scan_patt = sequence_points(xo, pattern_type)
            num_pts = numel(xo);
            num_rows = size(xo,2);
            switch pattern_type
                case 'serial_y_scan'
                    scan_patt = 1:(num_pts);
                % case 'serial_x_scan'            
                %     elem_num = reshape(1:(num_pts), size(xo));
                %     scan_patt = zeros(1, num_pts);
                %     for row_i = 1:num_rows
                %         % kept this bit of code here for historical reasons
                %         % (to see what had been done in the past), but
                %         % somewhere along the line the setting of patt_size
                %         % has been lost...
                %         s_i = (patt_size*(row_i - 1))+1;
                %         e_i = (patt_size*row_i);
                %         scan_patt(s_i:e_i) = elem_num(row_i,:);
                %     end
                case 'meander_x_scan'
                    elem_num = reshape(1:(num_pts), size(xo));
                    scan_patt = zeros(1, num_pts);
                    odd_row = true;
                    for row_i = 1:num_rows 
                        s_i = (num_rows*(row_i - 1))+1;
                        e_i = (num_rows*row_i);
                        if odd_row
                            scan_patt(s_i:e_i) = elem_num(row_i,:);
                        else
                            scan_patt(s_i:e_i) = elem_num(row_i,end:-1:1);
                        end
                        odd_row = not(odd_row);
                    end
                otherwise
                    error('unknown scan type');
            end        
        end
        
        function [probe_cmp_mag, probe_phs, out] = perfect_probe(dpSize, recon_pixel_size, ...
            defocus, beam_half_angle, lambda, samp_factor, Cs)
        
            % Here we follow the sign convention where negative defocus
            % produces and underfocused probe, which in the case of a
            % conventional STEM geometry puts a cross over below the sample
            % plane, and overfocus (+Ve) defocus produces a cross over
            % above.
            %
            % As a simple check of sign conventions here, if there is
            % cross over above the sample (+Ve defocus) then the diffraction 
            % pattern should look like a shadow image of the object.
            % 
            % Example usage:
            %     PS = ptychoset_mult;
            %     test_prb = PS.perfect_probe(512, 25e-12, 200e-9, 0.013, 2.5e-12, 1, 0);
            %     test_prb(:,1:256) = 0.5*test_prb(:,1:256); 
            %     make left side darker than right.
            %     figure; imagesc(abs(test_prb));
            %     figure; imagesc(abs(ft2(test_prb,1)));
            %
            % Should produce a diffraction pattern where left side is
            % darker than the right.
            %
            % Also the convention here, is set up so that negative defocus
            % can 'compensate' for +Ve spherical aberration as usual.
            
            
            if nargin < 7 || isempty(Cs)
                Cs = 0;
            end
        
            % quick hack trial of including Cs for recent recipes. This
            % needs to be integrated properly into the recipes/ experiments
            % of the driving spreadsheet etc.
            ps = @(alpha) (2*pi/lambda) * (((defocus/2) * alpha.^2) + ((1/4)*Cs*alpha.^4));
%             figure; atemp = linspace(0, beam_half_angle, 100); plot(atemp, ps(atemp)); 
%             ps = @(alpha) (2*pi/lambda) * ((defocus/2) * alpha.^2);

            if nargin < 6 || isempty(samp_factor)
                samp_factor = 1; 
            end

            [x, y] = meshgrid(-(dpSize/2):((dpSize/2)-1));
            r = sqrt(x.^2 + y.^2);
            r_m = r*recon_pixel_size;

            r_step = samp_factor*recon_pixel_size; % used to speed things 
            % up a bit. Interpolation should be fine in most situations.

            probe_image_fov = dpSize * recon_pixel_size;
            r_vals = 0:r_step:((probe_image_fov / 2)*sqrt(2));
            psi = zeros(1,length(r_vals));
            % mebs type approach..
            %{
            integrand_R = @(a, r) exp(-1i*ps(a)) .* besselj(0, (2*pi/lambda) .* r .* a) .* a;
            for ii = 1:length(r_vals)
                
                int_A = @(a) integrand_R(a, r_vals(ii));
                
                psi(ii) =  integral(int_A, 0, beam_half_angle,...   
                           'RelTol',1e-8,'AbsTol',1e-13);
            end
            psi = (2*pi/lambda)*sqrt(1/(pi*beam_half_angle^2))*psi;
%             psi = conj(psi);
            
            %}
            %{r
            integrand_R = @(r, rho) r .* exp(-1i*ps(r)) .* besselj(0, (2*pi .* rho .* r));
            for ii = 1:length(r_vals)
                
                int_A = @(r) integrand_R(r, r_vals(ii)/lambda);
                
                psi(ii) =  exp(1i*(pi/lambda)*r_vals(ii).^2) * ...
                           integral(int_A, 0, beam_half_angle,...   
                           'RelTol',1e-8,'AbsTol',1e-13);
            end
            % additional phase factor is usually pretty tiny, but might
            % come in at somepoint with low energy beam?...
            % Above eqn is also just rewritten to make it easier to tie up
            % with standard texts, such as goodman's fourier optics.
            psi = (2*pi/lambda)*sqrt(1/(pi*beam_half_angle^2))*psi;
            %}
            % in the early days of this code, when I had no
            % aberrations other than defocus to be concerned about, I 
            % inverted the sign of the angle by changing the sign of the
            % defocus
            
            if nargout == 3
                % just give radial profile, and not the 2D map.
                % Also do not normalize (to allow creation beam profiles).
                out.r = r_vals;
                out.psi = psi;
                % also compute the radial contained current integral:
                out.current = cumtrapz(2*pi*r_vals.*abs(out.psi).^2);                
                probe_cmp_mag = [];
                probe_phs = [];
            elseif nargout < 3
                psi = psi/abs(psi(2)); % normalize
                probe_cmp_mag = reshape(interp1(r_vals, abs(psi), r_m(:), 'nearest'), size(r));
                probe_phs = reshape(interp1(r_vals, angle(psi), r_m(:), 'nearest'), size(r));
            end
            
            if nargout == 1
                probe_cmp_mag = probe_cmp_mag.*exp(1i*probe_phs);
                probe_cmp_mag(isnan(probe_cmp_mag)) = 0;
            end
            

        end    

        function [probe_modulus, src_struct] = probe_guess_from_dm3(...
                source_image, src_type, ...
                output_pixel_size, output_num_pix, rot_angle_deg  )


            % output_pixel_size = 3.9061e-11;
            % output_num_pix = 1024;  % must be even!
            % rot_angle_deg = 0;

            % from ptycho note:
            % recon_pixel_size = lambda/(angular_range);

            if ischar(source_image)
                switch src_type
                    case 'dm3'
                        % prb_im = importDM(source_image); From Mike
                        % im = prb_im.TheData.Data;
                        % [per_pix, ~] = prb_im.TheData.Calibrations.Dimensions.Scale;
                        % [im_unit, ~] = prb_im.TheData.Calibrations.Dimensions.Units;
                        prb_im = DM3Import( source_image ); % From Robb
                        im = prb_im.image_data;
                        per_pix = prb_im.xaxis.scale;
                        im_unit = prb_im.xaxis.unuts;
                    case 'mat'
                        prb_im = load(source_image);
                        im = prb_im.Image;
                        per_pix = prb_im.Scale;
                        im_unit = prb_im.Units;
                    otherwise
                        error('Unknown source file type');
                end
            else
                if isstruct(source_image)
                    im = source_image.Image;        
                    per_pix = source_image.Scale;
                    im_unit = source_image.Units;
                else
                    error('Unrecognized variable passed');
                end
            end

            probe_image_fov = output_num_pix * output_pixel_size;


            switch im_unit
                case 'nm'
                    multip = 1e-9;
                case 'ï¿½m'
                    multip = 1e-6;
                case 'um'
                    multip = 1e-6;
                case 'micron'
                    multip = 1e-6;
                otherwise
                    error('Unknown unit');
            end

            % Recentre the probe image:
            cent_pos = FindCentre(im);
            im = abs(Recentre(im,cent_pos));

            raw_pix_size = per_pix * multip;
            num_of_pix_in_fov = probe_image_fov/raw_pix_size;
            scale_factor = output_num_pix/num_of_pix_in_fov;

            % crop image so we work only part of it that is needed:
            cropped_pix_size = 2*ceil((num_of_pix_in_fov * 1.5)/2);
            base_size = size(im);

            r_rng = ((base_size(1)-cropped_pix_size)/2 + 1):...
                   ((base_size(1)-cropped_pix_size)/2 + cropped_pix_size) ;
            c_rng = ((base_size(2)-cropped_pix_size)/2 + 1):...
                   ((base_size(2)-cropped_pix_size)/2 + cropped_pix_size) ;    

            cropped_raw = im(r_rng, c_rng);

            % scale this
            scaled_im = imresize(cropped_raw, scale_factor);

            isodd = @(x) mod(x,2);

            if isodd(size(scaled_im,1))
                scaled_im = scaled_im(1:(end-1),1:(end-1));
                scaled_im = abs(Recentre(scaled_im,FindCentre(scaled_im)));
            end


            % rotate it:
            rot_im = imrotate(scaled_im, rot_angle_deg, 'bilinear','crop');

            % final crop
            base_size = size(rot_im);

            r_rng = ((base_size(1)-output_num_pix)/2 + 1):...
                   ((base_size(1)-output_num_pix)/2 + output_num_pix) ;
            c_rng = ((base_size(2)-output_num_pix)/2 + 1):...
                   ((base_size(2)-output_num_pix)/2 + output_num_pix) ;  

            rot_im = scaled_im(r_rng, c_rng);
            probe_modulus = sqrt(rot_im);
            % figure; imagesc(probe_modulus);axis image; colormap jet;

            if nargout == 2
                src_struct = struct('Image', im , 'Scale', per_pix, 'Units', im_unit);
            end

        end
        
        function [drift_mag_nm, drift_ang_deg, drift_mag_nmpersec] = calc_drift(first_image, last_image)

            % first_image = ['H:\Arthur\201611xx_Ptycho\',data_setname,'\camera_setup\before exposure.dm3'];
            % last_image = ['H:\Arthur\201611xx_Ptycho\',data_setname,'\camera_setup\after exposure.dm3'];

            im1 = DM3Import(first_image);
            im1_Data = im1.image_data;
            [~,~,EXT] = fileparts(last_image);
            if strcmpi(EXT,'.tif')
                t = Tiff(last_image,'r');
                t_desc = getTag(t,'ImageDescription');
                k = strfind(t_desc,'Time=');
                timeread = t_desc((k+5):(k+5+7));
                timestr = 'HH:MM:SS';
                after_time = datenum(timeread,timestr);
                im2_Data = t.read;
                t.close;
                tif_read = true;
            else
                im2 = DM3Import(last_image);
                im2_Data = im2.image_data;
                tif_read = false;
            end
            units_per_pix = im1.xaxis.scale;
            im_unit = im1.xaxis.units;
     
            switch im_unit
                case 'nm'
                    multip = 1;
                case 'ï¿½m'
                    multip = 1000;
                case 'um'
                    multip = 1000;
                case 'micron'
                    multip = 1000;
                otherwise
                    error('Unknown unit');
            end

            nm_per_pix = units_per_pix * multip;
            
            sub_pix = 10;
            a = dftregistration(fft2(im2_Data) , fft2(im1_Data) , sub_pix); 

            % df_px.x = -a(4);
            % df_px.y = a(3);

            drift_vec = complex(-a(4), a(3))*nm_per_pix;
            drift_mag_nm = abs(drift_vec);
            drift_ang_deg = (360/(2*pi))*angle(drift_vec);
            if nargout == 3

                dm3datestr ='ddd mmm dd HH:MM:SS yyyy';
                d1 = datevec(datenum(char(im1.AcquisitionStartTimeStamp)',dm3datestr));
                if tif_read
                    d2 = datevec(after_time);
                    % ignore the yr, day, and month as currently that's
                    % not in meta info
                    d1(1:3) = [2020,1,1];
                    d2(1:3) = [2020,1,1];
                else
                    d2 = datevec(datenum(char(im2.AcquisitionStartTimeStamp)',dm3datestr));
                end
                
                elapsed_time = etime(d2, d1);
                drift_mag_nmpersec = drift_mag_nm/elapsed_time;
            end
        end
        
        function out_func = create_sfits_from_coefs(fitdat)

           for ii = 1:length(fitdat.fnames)
               c.(fitdat.fnames{ii}) = fitdat.vals(ii);
           end
           
           switch fitdat.type
               case 'poly22'
                   out_func = @(x,y) c.p00 + c.p10*x + c.p01*y + c.p20*x^2 + c.p11*x*y + c.p02*y^2;
               case 'poly11'
                   out_func = @(x,y) c.p00 + c.p10*x + c.p01*y;
               otherwise
                   disp('unknown fit type');
           end
            
        end
        
        function [out_func, fitdat] = conv_sfits_to_poly(sf)
           fnames = coeffnames(sf);
           vals = coeffvalues(sf);
           for ii = 1:length(fnames)
               c.(fnames{ii}) = vals(ii);
           end
           switch type(sf)
               case 'poly22'
                   out_func = @(x,y) c.p00 + c.p10*x + c.p01*y + c.p20*x^2 + c.p11*x*y + c.p02*y^2;
               case 'poly11'
                   out_func = @(x,y) c.p00 + c.p10*x + c.p01*y;
               otherwise
                   disp('unknown fit type');
           end
           fitdat.fnames = fnames;
           fitdat.vals = vals;
           fitdat.type = type(sf);
    
        end 
                 
        function p = parse_and_check_ip(param_in, input_info)
        % parses, checks the input and sets to default where needed.
        param_names = input_info(:,1);
        param_defaults = input_info(:,2);

        % If there happens to be addional field on param_in, let them through
        % unparsed. This is to deal with situation where a struct is passed
        % directly in.
        p = param_in;
        % check all the ones that are listed.

        for p_ind = 1:length(param_names)
           if isfield(param_in,param_names{p_ind})
               arg_in = param_in.(param_names{p_ind});
               if ischar(arg_in)
                   % Then the input has come either from the command line or is a
                   % text input
                   arg = str2num(arg_in);  %#ok<ST2NM>
                   % if the characters represent a number or an array this will
                   % convert, otherwise it will be empty, in which case the
                   % argument should be copied as is - but stripped of any possible
                   % leading or trailing whitespace:
                   if isempty(arg)
                       arg = sscanf(arg_in,'%s');
                   end
                   % now the problem of how to represent empty or default elements
                   % as text input
                   if ischar(arg) && ...
                           (strcmpi(arg,'empty') || ...
                            strcmpi(arg,'none') )
                   arg = [];
                   end 
                   if ischar(arg) && ...                    
                            (strcmpi(arg,'default') || ...
                             isempty(arg))
                   arg = param_defaults{p_ind};
                   end
               else
                   arg = arg_in;
               end
               p.(param_names{p_ind}) = arg;
           else
               p.(param_names{p_ind}) = param_defaults{p_ind};
           end 
        end
    end
        
    
    end

end

function [ims, tags] = read_dm3stack(filenames)

    num_files = length(filenames);
    ims = cell(1,length(num_files));
    tags = cell(1,length(num_files));
    for ii = 1: num_files
        imname = filenames{ii};
        read_data = DM3Import(imname);
        ims{ii} = single(read_data.image_data);
        tags{ii} = rmfield(read_data, 'image_data');             
    end
end

