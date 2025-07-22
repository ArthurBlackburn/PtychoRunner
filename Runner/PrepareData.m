function PrepareData(varargin)
% A function for preparing diffraction data from hpl files produced from
% Azorus. Performs recentering, and reordering if necessary. No doubt this
% could be made run alot quicker with python, but for now this will do the
% job...
%
% Examples
%
% See "TestScripts\Part_0_Prepare_Data.m" for other examples of usage.
%
% PrepareData('input_file_or_dir','D:\Data\Pt_TEST\AuTest.hp');
%
% PrepareData('input_file_or_dir','D:\Data\Pt_TEST\AuTest.hp',...
% 'positions_trans',struct('trans_coords',5,'rot_angle_deg',0),'centre_pos',[220, 287]);
%
% PrepareData('input_file_or_dir','D:\Data\Arthur\20230622\AuTest.hp',...
% 'positions_trans',struct('trans_coords',5,'rot_angle_deg',0),'centre_pos',[220, 287],...
% 'outputDPSizes',{256,512},'determine_masks',false);
% 
% PrepareData('input_file_or_dir','D:\Data\Arthur\20230622\PtIr01.hp',...
% 'positions_trans',struct('trans_coords',5,'rot_angle_deg',0),...
% 'outputDPSizes',{256,512,1024},'determine_masks',false);
%     >> 
%     Accept x: 213, y: 289 as position (y/n)?
%     y
%     Centre of disc taken as (Row #, Column #): (289, 213) [green x]
%     i.e. 'centre_pos',[213, 289]
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

import plotting.*   % from ptychoshelves, for imagesc3D.

FileOrDir = @(x) (isfile(x) || isfolder(x));

in = inputParser;
in.addParameter('do_all_in_dir',false, @(x) assert(islogical(x)));
in.addParameter('input_file_or_dir', [], @(x) assert(FileOrDir(x),'Input file or directory must be specified'));
in.addParameter('save_dir',[], @(x) assert(~isempty(x),'Output Directory must be specified'));
in.addParameter('interactive',true);
in.addParameter('centre_pos',[]);
% in.addParameter('centre_pos',[col_0, row_0]);
in.addParameter('determine_masks',true);
in.addParameter('do_plots',true);
in.addParameter('disp_gamma',0.4);
in.addParameter('example_num',[]);
in.addParameter('default_centre',[]);
in.addParameter('positions_trans',[]); % makes user input values
% in.addParameter('positions_trans',struct('trans_coords',5,'rot_angle_deg',0));
in.addParameter('outputDPSizes',{'same'});
% in.addParameter('outputDPSizes',{'same',256,1024});
in.addParameter('pattern_size_out',[]);
% in.addParameter('pattern_size_out',[40,40]);

nm_V = 1; % leave this at one, so that nm_V can be set in excel sheet / elsewhere


in.parse(varargin{:});
p = in.Results;

if ~isempty(p.save_dir)
    try
        if ~isfolder(save_dir)
            mkdir(p.save_dir);
        end
    catch
        warning('Problem creating save directory. Placing output files in same directory as input');
        p.save_dir = [];
    end
end

if isempty(p.save_dir)
    [dirname, ~, ext] = fileparts(p.input_file_or_dir);
    if isempty(ext) && isfolder(p.input_file_or_dir)
        p.save_dir = p.input_file_or_dir;
    else
        p.save_dir = dirname;
    end
    if isfolder(p.save_dir)
        warning('Save directory was not specified. Will place output in same directory as input.');        
    end
end

%% init filenames
if p.do_all_in_dir
    dir_fnames = dir(fullfile(p.input_file_or_dir,'*.hp'));
    dir_fnames = struct2cell(dir_fnames);
    all_hps = dir_fnames(1,:);
else
    all_hps = {p.input_file_or_dir};
end

%%

first_run = true; % usually this is set true so the the first run is interactive.
% first_run = false;


%%
for ii = 1:length(all_hps)

    fname = all_hps{ii};
    dat_in = loadHPL(fname);
    numDPs = size(dat_in.diffraction.micrograph,1);
    sourceDPsize = size(dat_in.diffraction.micrograph,[2 3]);

    if prod(dat_in.attrs.meshParams.shape) ~= numDPs
        % then likely this is data acquired prior to early 2021, 
        % which had Robb's error in the pattern shape
        pattern_shape = dat_in.attrs.meshParams.shape.' + 1; %
        fprintf('Corrected Pattern shape by adding 1 to dimensions');
    else
        pattern_shape = dat_in.attrs.meshParams.shape.';
    end
    
    fprintf('Pattern Shape: %d x %d\n',pattern_shape);
    
    if prod(pattern_shape) ~= numDPs
        warnmsg = ['For rectangular scans, product of pattern shape should agree with number of DPs.',...
            'This is not the case here, so problems / errors might arise later if trying to make a subset.'];
        warning(warnmsg)
    end
    
    while ~isempty(p.pattern_size_out) && ~all(p.pattern_size_out <= pattern_shape)
        fprintf(['pattern size out: %d, %d is inconsistent with\n',...
                 'pattern size in: %d, %d\n'],...
                 p.pattern_size_out, pattern_shape);
        try
            p.pattern_size_out = input('Enter output pattern shape as array: ');
        catch ME
            warning(ME.message);
            fprintf('Retry entry\n');
        end             
    end

    if isempty(p.example_num)
        p.example_num = round(rand(1)*numDPs);
    end
    if isempty(p.default_centre)        
        p.default_centre = floor(size(dat_in.diffraction.micrograph,[2 3])/2)+1;
    end

    % ---------------------------------------------------------------------
    
    for osi = 1:length(p.outputDPSizes)
        outDPsize = p.outputDPSizes{osi};
                
        if ischar(outDPsize) && strcmpi(outDPsize,'same')
            newDPSize = sourceDPsize;
            noresize = true;
        else
            if outDPsize == sourceDPsize
                noresize = true;
            else
                noresize = false;    
            end
            newDPSize = outDPsize;            
        end
             
        voltage_coords = dat_in.coords;

        if first_run || p.interactive
            if isempty(p.positions_trans)
                p.positions_trans = get_posn_transforms(voltage_coords, nm_V, fname);
            end
            if isempty(p.centre_pos)
                % Calc C of M:
                p.centre_of_mass = FindCentre(squeeze(sum(dat_in.diffraction.micrograph,1)) );
                display_dps();
                p.centre_pos  = get_DP_centre();
            end
            first_run = false;
        end
        
        [dp, data_region] = recentre_and_resize(newDPSize);
        
        if p.do_plots
            p.default_centre = floor(size(dp,[2 3])/2)+1;
            display_dps(dp, sprintf('recentred, %s, [%d, %d] px', strrep(fname,'\','\\'), newDPSize, newDPSize));
            hold on;
            plot([data_region(3), data_region(4), data_region(4), data_region(3), data_region(3)],...
                 [data_region(2), data_region(2), data_region(1), data_region(1), data_region(2)], 'LineWidth',2, 'LineStyle','--');
        end

        if p.determine_masks
            outliers = get_outlier_pixel_masks();
            fona = fieldnames(outliers);
            for ona = 1:numel(fieldnames(outliers))
                oname = fona{ona};
                if ismatrix(outliers.(oname)) && all(size(outliers.(oname)) == sourceDPsize)
                    o_in(1,:,:) = outliers.(oname);
                    o_out = recentre_and_resize(newDPSize, o_in);
                    outliers.(oname) = squeeze(o_out);                 
                end
            end           
        end
        
        dp = permute(dp,[2 3 1]); % put in the order required for other programs...
         
        % set up the names to save outputs to:
        [~,fn,~] = fileparts(fname);
        savenames.suffix = sprintf('_rc_%d',newDPSize(1)); % '_rc_' is meant to signify that the the data has been recentred.
        savenames.matlab_shifted_and_cropped  = fullfile(p.save_dir, [fn,savenames.suffix,'.mat']);
        savenames.coords = fullfile(p.save_dir, [fn,'_COORDS.mat']);
        
        if p.determine_masks
            savenames.matlab_shifted_and_cropped_mask = fullfile(p.save_dir, [fn, savenames.suffix,'_MASK.mat']);            
        else
            savenames.matlab_shifted_and_cropped_mask = [];
        end
        
        if isempty(p.pattern_size_out) || all(p.pattern_size_out == pattern_shape)
            savenames.dp_datafile_name = 'none';
            savenames.out_posfilename = fullfile(p.save_dir, [fn,savenames.suffix,'_position.hdf5']);
            savenames.out_dpfilename = 'none';
        else
            savenames.dp_datafile_name = savenames.matlab_shifted_and_cropped;
            dset_basename = sprintf('%s_%dby%d_%d',fn, p.pattern_size_out(1), p.pattern_size_out(2), newDPSize(1));
            savenames.out_posfilename = fullfile(p.save_dir, [dset_basename,'_position.hdf5']);
            savenames.out_dpfilename = fullfile(p.save_dir, [dset_basename,'_dp.mat']);            
        end
                        
        % set up vars for save
        centre_pos = p.centre_pos;
        positions_trans = p.positions_trans;
        % save the dp stack
        save(savenames.matlab_shifted_and_cropped,'dp','data_region','pattern_shape','positions_trans','centre_pos','voltage_coords','-v7.3');
        % also save the coords in separate file at the moment.
        save(savenames.coords , 'voltage_coords', 'pattern_shape', 'positions_trans');

        if p.determine_masks
            save(savenames.matlab_shifted_and_cropped_mask,'outliers','data_region','pattern_shape','-v7.3');
        end

        clearvars dp

        create_ptyshelves_subset( savenames.dp_datafile_name,...
                                  savenames.coords,...
                                  savenames.out_dpfilename,...
                                  savenames.out_posfilename,...
                                  pattern_shape, p.pattern_size_out,...
                                  nm_V,  p.positions_trans.rot_angle_deg, p.positions_trans.trans_coords,...
                                  savenames.matlab_shifted_and_cropped_mask, false)
    end                          

end

% ---- NESTED FUNCTIONS -----
    function display_dps(plotdp, title_suffix)
        
        if nargin < 1
            plotdp = dat_in.diffraction.micrograph;
        end
        if nargin < 2
            title_suffix = [];
        end
        
        % display diffraction patterns:    
        figure; 
        imagesc3D(permute(plotdp,[2 3 1]));
        ts = ['DP Stack ',title_suffix];
        axis image
        colorbar;
        title(ts,'Interpreter','none');
        hold on;

        % overlay the default centre and the centre of mass of the stack    
%         plot(p.default_centre(2),p.default_centre(1),'bx');
        plot(p.default_centre(2),p.default_centre(1),'go');
        legtext = {'default centre'};

        centre_of_m = FindCentre(squeeze(sum(plotdp,1)) );
%         plot(centre_of_m(2), centre_of_m(1),'mx');
        plot(centre_of_m(2), centre_of_m(1),'mo');
        legtext{end+1} = 'centre of mass';    
        legend(legtext);% 'Position',

        drawnow;                
    end


    function centre_pos  = get_DP_centre()            
        accept_res = 'n';

        centre_select_method = input(sprintf([...
                'Manually select centre position of DP (m), \n',...
                'use default stored value [R:%d, C:%d] [shown blue] (d) ,\n',...
                'or use centre of mass of stack to shift all [shown pink] (c)',...
                '?\n(Zoom plot before continuing). (m/d/c)>'],...
                 p.default_centre(1), p.default_centre(2)),'s');


        while ~strcmpi(accept_res,'y')
            bad_entry = false;
            switch centre_select_method
                    case 'm'
                        fprintf('Select centre position\n:');
                        centre_pos = round(ginput(1));
                        centre_pos = centre_pos([2 1]);           
                    case 'd'
                        centre_pos = p.default_centre;
                    case 'c'
                        centre_pos = centre_of_mass;
                    otherwise
                        fprintf('Incorrect entry. Try again!\n');
                        bad_entry = true;
            end     

            if ~bad_entry
                c_marker = plot(centre_pos(2), centre_pos(1), 'kx', 'MarkerSize', 20); %#ok<NASGU>
                drawnow;
                accept_res = input(sprintf('Accept x: %d, y: %d as position (y/n)?\n',...
                    centre_pos(1), centre_pos(2)),'s');
                if strcmpi(accept_res, 'y')
                    % Now carry on get radius:
                    fprintf('Centre of disc taken as (Row #, Column #): (%d, %d) [green x]\n',...
                    centre_pos(2), centre_pos(1));   
                else
                    accept_res = 'n';
                    delete c_marker
                end
            end
        end
    end

    function outliers = get_outlier_pixel_masks()

        exampleDP = squeeze(dat_in.diffraction.micrograph(p.example_num,:,:));

        % See if any pixels are above the expected max count:
        hot_count = 2^16-1;
        outliers.hot_pix_examp = exampleDP >= hot_count;                
        outliers.sum_hot_mask = squeeze(sum(dat_in.diffraction.micrograph >= hot_count, 1));
        outliers.hot_some_of_the_time = outliers.sum_hot_mask > 0;

        hot_thresh = graythresh(outliers.sum_hot_mask)*max(outliers.sum_hot_mask(:));
        outliers.thresh_hot_mask = outliers.sum_hot_mask > hot_thresh;


        % See if any pixels are always zero:               
        cold_count = 0;
        outliers.cold_pix_examp = exampleDP == cold_count;  
        outliers.sum_cold_mask = squeeze(sum(dat_in.diffraction.micrograph == cold_count, 1));                                
        outliers.cold_some_of_the_time = outliers.sum_cold_mask > 0;


        % there are quite a lot that are not cold all the time: these must be ones
        % outside aperture that occasionally receive stray electrons.
        % lets make a threshold using Otsu method:
        cold_thresh = graythresh(outliers.sum_cold_mask)*max(outliers.sum_cold_mask(:));
        outliers.thresh_cold_mask = outliers.sum_cold_mask > cold_thresh;

        outliers.suggested_bad_pixels = outliers.thresh_hot_mask | outliers.thresh_cold_mask;

        % figure; 
        % #2
        if p.do_plots
            figure;
            % #1
            ax11 = subplot(2,3,1);
            imagesc(exampleDP.^p.disp_gamma); axis image;
            title(sprintf('Example DP #%d, Gamma = %3.3f', p.example_num, p.disp_gamma));
            colorbar;

            ax12 = subplot(2,3,2);
            imagesc(outliers.suggested_bad_pixels); axis image; colorbar
            title('Suggested Bad Pixels');

            ax13 = subplot(2,3,3);
            imagesc(outliers.sum_hot_mask/size(dat_in.diffraction.micrograph,1)); axis image; colorbar
            title('Fraction of frames where pixel appears hot.');

            ax14 = subplot(2,3,4);
            imagesc(outliers.sum_cold_mask/size(dat_in.diffraction.micrograph,1)); axis image; colorbar
            title('Fraction of frames where the pixel is cold');

            ax15 = subplot(2,3,5);
            imagesc(outliers.thresh_hot_mask); axis image; colorbar;
            title('Mask used to mark hot pixels that c/should possibly be ignored');

            ax16 = subplot(2,3,6);
            imagesc(outliers.thresh_cold_mask); axis image; colorbar;
            title('Mask used to mark cold pixels that c/should possibly be ignored');

            linkaxes([ax11 ax12 ax13 ax14 ax15 ax16]);
            sgtitle(sprintf('Plots Related to %d cropped data',newDPSize));
        end
    end

    function [imout, data_region_out] = crop_pad_hpl(Nout, imgin, data_region_in)
        outside_value = 0;
        
        if nargin < 2
            imgin = dat_in.diffraction.micrograph;
        end
        Nframes = size(imgin,1);
        Nin = size(imgin,[2 3]);
        
        if nargin < 3
            data_region_in = [1, Nin(1), 1, Nin(2)];
        end
        
        if numel(Nout) == 1
            Nout = Nout*[1 1];
        end
                
        center = floor(Nin(1:2)/2)+1;

        imout = zeros([Nframes, Nout],'like',imgin) + outside_value;
        centerout = floor(Nout/2)+1;

        cenout_cen = centerout - center;
        
        R1_o = max(cenout_cen(1)+1,1);
        R2_o = min(cenout_cen(1)+Nin(1),Nout(1));
        C1_o = max(cenout_cen(2)+1,1);
        C2_o = min(cenout_cen(2)+Nin(2),Nout(2));
        data_region_o = [R1_o, R2_o, C1_o, C2_o];
        
        
        R1_i = max(-cenout_cen(1)+1,1);
        R2_i = min(-cenout_cen(1)+Nout(1),Nin(1));
        C1_i = max(-cenout_cen(2)+1,1);
        C2_i = min(-cenout_cen(2)+Nout(2),Nin(2));
        
        data_region_i = [R1_i, R2_i, C1_i, C2_i];
        
        imout(:,...
              R1_o : R2_o,...
              C1_o : C2_o ) = ...
        imgin(:,...
              R1_i : R2_i,...
              C1_i : C2_i );
          
        used_data_region = [max(R1_i,data_region_in(1)),...
                            min(R2_i,data_region_in(2)),...
                            max(C1_i,data_region_in(3)),...
                            min(C2_i,data_region_in(4))] - ...
                          data_region_i;
        data_region_out = data_region_o + used_data_region;

    end

    function [dpout, data_region] = recentre_and_resize(newsize, dpin, nosz)
        outside_value = 0;
        
        if nargin < 3
            nosz = noresize;
        end
        if nargin < 2
            dpin = dat_in.diffraction.micrograph;
        end
        if nargin < 1
            newsize = newDPSize;
        end        
        
        if nosz
            [dpout, ~, ~, data_region]  = recentre_HPL_dpstack_simple(dpin, p.centre_pos, outside_value);
        else            
            if all(newsize < sourceDPsize)
                [dpout, ~, ~, data_reg01]  = recentre_HPL_dpstack_simple(dpin, p.centre_pos, outside_value);                
                [dpout, data_region] = crop_pad_hpl(newsize, dpout, data_reg01);
            end
            if all(newsize > sourceDPsize)
                [dpout, data_reg01] = crop_pad_hpl(newsize, dpin);
                offset = floor(newsize/2) - floor(sourceDPsize/2); 
                new_centre = p.centre_pos + offset;
                [dpout,~,offset] = recentre_HPL_dpstack_simple(dpout, new_centre, outside_value);
                data_region = data_reg01+[offset(1), offset(1), offset(2), offset(2)];
            end            
        end
    end

end

function positions_trans = get_posn_transforms(voltage_coords, nm_V, fname)
    convert_to_angstrom = false; % also leave this alone
    
    accept_res = 'n';
    while ~strcmpi(accept_res,'y')
        figure;
        ax01 = subplot(1,2,1);
        plot(voltage_coords(:,1),voltage_coords(:,2));
        axis equal;
        xlabel('Voltage Coords, Column 1');
        ylabel('Voltage Coords, Column 2');
        xlim([-6 6]); ylim([-6 6]);
        hold on;
        plot(voltage_coords(1:10,1),voltage_coords(1:10,2),'ro');
        sgtitle(sprintf('Data File: %s',fname),'interpreter','none');
        title('Source Voltage Coords','interpreter','none');
        positions_trans.trans_coords = ...
        input(['Enter which transformation makes the voltage coordinates cartesian ( 1 - 8)\n',...
            '[i.e. which transformation puts first points (marked red) in top left for a _regular scan_\n',...
              'noting that the Azorus STEM display might by accident be transposed etc, from display on microscope GUI...]\n',...
                    'case 1, positions_real(x,y) = [-col 2,  col 1];\n',...
                    'case 2, positions_real(x,y) = [ col 2,  col 1];\n',...
                    'case 3, positions_real(x,y) = [ col 1,  col 2];\n',...
                    'case 4, positions_real(x,y) = [ col 1, -col 2];\n',...
                    'case 5, positions_real(x,y) = [ col 2, -col 1];\n',...
                    'case 6, positions_real(x,y) = [-col 2, -col 1];\n',...
                    'case 7, positions_real(x,y) = [-col 1,  col 2];\n',...
                    'case 8, positions_real(x,y) = [-col 1, -col 2];\n',...
                    '(1 - 8) > ']);
        positions_trans.rot_angle_deg = input('Angle that was applied in Azorus transformation of coordinates (degrees):\n');

        probe_positions_0 = gen_coords_for_ptychoshelves(voltage_coords, [],...
            nm_V, positions_trans.rot_angle_deg, convert_to_angstrom, ....
            positions_trans.trans_coords);
        ax02 = subplot(1,2,2);
        plot(probe_positions_0(:,1), probe_positions_0(:,2));
        xlim([-6 6]); ylim([-6 6]);
        axis equal;
        hold on;
        plot(probe_positions_0(1:10,1), probe_positions_0(1:10,2),'ro');            
        xlabel('Cartesian Coords, X');
        ylabel('Cartesian Coords, Y');
        title(sprintf('Coords after trans case %d',positions_trans.trans_coords),'interpreter','none');
        linkaxes([ax01 ax02]);
        accept_res = input('Is this transformation correct? (y/n)\n','s');
    end

end



