function [res, pairs] =  CalcFRCMeas(pairs, varargin)
% Wrapper function for calculating FRC measures calling on functions of PtychoShelves but
% using output obtained from PtychoRunner. It can calculate FRC over subregions as used within
% particulate_FRC_res, i.e. it calcs of FRC measures over a series of regions of interest in 
% over the field of view.
%
% Examples:
%
% r0 = load('D:\Data\Ptycho20kV_QuickDL\AuonC\S_CR_26\Au20kV_2022_B_CR_26_recipe015.mat');
% r1 = load('D:\Data\Ptycho20kV_QuickDL\AuonC\S_CR_26\Au20kV_2022_C_CR_26_recipe015.mat');
% 
% res = CalcFRCMeas({r0, r1});
%
% Gives output that starts:
% 
%         Manually select the compared region ... 
% >> (User selects ROI)
%         ===========================
%         Selected region: {935:3321,952:3356}
%         ===========================
%         Zoom figures now prior to picking features. Press any key when done
% >> (Do zooming manually on regions that look similar)
%         Click on a feature on figure 21
%         Click on a feature on figure 22
%         Initial Offset Guess: [x: 4, y: 12]
%
% Then the we can automatically search a region for best alignment, noting that the above gives 
% the compared region position 
% 
% {Row Lowest Num: Row Highest Num, Column Num Left-Most : Column Num Right-Most}
%
% which should be fed into the FOV_corners param here, where the order of points is
%
% [Row Lowest  Num, Column Num Left-Most , Row Highest Num, Column Nun Right-Most]
%
% e.g. for the above 
% params = {'passtoalign',{'FOV_corners',[935,952,3321,3356],'guessx',4,'guessy',12,...
%           'searchForBest',true, 'searchGridPts',5, 'autoplotclose', false}};
%
% res = CalcFRCMeas({r0, r1}, params{:});
%
% See further examples in particulate_FRC_res.m
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

SetupPathsN(true);
import utils.*

in = inputParser;
in.addParameter('pair_directories',[]);
in.addParameter('s_modifiers',[]);
in.addParameter('saveresname',[]);
in.addParameter('passtoalign',{});
in.parse(varargin{:});
p = in.Results;

if ~isempty(p.saveresname)
    p.passtoalign{end+1} =  'saveresname';
    p.passtoalign{end+1} = p.saveresname;
end

if isempty(pairs)

    set0_to_set1_modifier = @(name0, pn) strrep(name0, p.s_modifiers{pn,1}, p.s_modifiers{pn,2});

    search_ext = '*.mat';

    pairs = cell(1,2);
    mp = 0;
    for pair_num = 1:size(p.pair_directories,1)

        set0 = dir(fullfile(p.pair_directories{pair_num, 1}, search_ext));
        set1 = dir(fullfile(p.pair_directories{pair_num, 2}, search_ext));
        name_mask = strcmpi('name',fieldnames(set1));
        set0 = struct2cell(set0);
        set0 = set0(name_mask,:);
        set1 = struct2cell(set1);
        set1 = set1(name_mask,:);

        num_set0 = length(set0);
        for ii = 1:num_set0
            if any(strcmpi(set0_to_set1_modifier(set0{ii},pair_num), set1))
                mp = mp+1;
                pairs{mp,1} = fullfile(p.pair_directories{pair_num, 1}, set0{ii});
                pairs{mp,2} = fullfile(p.pair_directories{pair_num, 2}, set0_to_set1_modifier(set0{ii}, pair_num));
            end
        end

    end
end

    
if size(pairs,1) > 1
    res = cell(size(pairs,1),1);
    for ii = 1:size(pairs,1)
        res{ii} = getFSCmeasure(pairs{ii,:}, p.passtoalign{:});
    end
else
    res = getFSCmeasure(pairs, p.passtoalign{:});
end

    
end

%%

function rrr = getFSCmeasure(pair, varargin)
import utils.*

    p = inputParser;
    p.addParameter('saveresname',[]);
    p.addParameter('savefigname',[]);
    p.addParameter('savefig',false);
    p.addParameter('manual_crop',true); %%
    p.addParameter('plot_level',4);
    p.addParameter('autoplotclose',false);
    p.addParameter('image_prop','complex');
    p.addParameter('do_align',true); % added to allow non-aligned use of subfields in particle_FRC_res
    p.addParameter('searchForBest',false);
    p.addParameter('roisize',[1024 1024]);
    p.addParameter('FOV_corners',[500,500,2000,2000]);
    p.addParameter('searchGridPts',6);
    p.addParameter('guessx',0); % guess x = round(xin0-xin1) where xin is a point(x) in image 0, corresponding to point(x) in image 1; (x = cols co-ordinate);
    p.addParameter('guessy',0); % guess y = round(yin0-yin1) where xin is a point(y) in image 0, corresponding to point(y) in image 1; (y = cols co-ordinate);
    p.addParameter('use_window',false,@islogical);
    p.addParameter('window_type','hann',@(x) ismember(x,{'hann','hamming'}));
    p.addParameter('mask_form',false);
    p.addParameter('crop','');
    p.addParameter('show_2D_fourier_corr',true);

    p.parse(varargin{:});
    rp = p.Results;
    
    [im0, im0_dat] = load_check_rotate(pair{1});
    [im1, im1_dat] = load_check_rotate(pair{2});
    if rp.mask_form
        im0 = make_and_mult_mask(im0);
        im1 = make_and_mult_mask(im1);
    end
    
    assert(im0_dat.recon_pixel_size == im1_dat.recon_pixel_size, 'Files must have same recon pixel size');
    assert(im0_dat.lambda == im1_dat.lambda, 'Files must have same wavelength');


    param.do_align = rp.do_align;
    param.lambda = im0_dat.lambda;
    param.pixel_size = im0_dat.recon_pixel_size;
    param.use_window = rp.use_window;
    param.window_type = rp.window_type;
    param.crop = rp.crop;

    if isempty(rp.savefigname) && rp.savefig
        rp.savefigname =  strrep(pair{1},'.mat','_complex_paper.fig');
    end

    param.plotting = rp.plot_level;
    % 'complex', 'phasor', 'phase', 'variation'
    % param.image_prop = 'phasor';
%     param.image_prop = 'complex';
    param.image_prop = rp.image_prop;
    param.prop_obj = false;
    if rp.searchForBest
        rp.manual_crop = false;
    end
%     if rp.manual_crop && ~isempty(param.crop)
    if rp.manual_crop
        param.crop = 'manual';
        param.GUIguess = true;        
    else
        param.crop = [];
        param.GUIguess = false;
        param.guessx = rp.guessx;
        param.guessy = rp.guessy;
    end
    param.verbose_level = 4;
    param.xlabel_type = 'resolution';
    param.SNRt = 0.2071;
    param.show_2D_fourier_corr = rp.show_2D_fourier_corr;

    if rp.manual_crop
        accept_res = 'n';
        while(strcmpi(accept_res,'n'))
            [rrr.resolution, rrr.stat, rrr.FSC, rrr.T, rrr.freq, rrr.n] = aligned_FSC(im0, im1, param);
            accept_res = input('Accept current alignment etc.? (y/n)\n','s');
            if strcmpi(accept_res,'n')
                fprintf('Starting again!\n');
                close('all');
            end
        end
    else
        rrr = struct('resolution', cell(rp.searchGridPts*[1 1]), 'stat', cell(rp.searchGridPts*[1 1]),...
            'mask',zeros(size(im0),'logical'));

        row_width = rp.FOV_corners(3) - rp.FOV_corners(1);
        col_width = rp.FOV_corners(4) - rp.FOV_corners(2);
        
        scan_size = [row_width, col_width] - rp.roisize;
        scan_step = scan_size/(rp.searchGridPts-1);
        Rstarts = round(rp.FOV_corners(1) + ((0:1:(rp.searchGridPts-1))*scan_step(1)));
        Cstarts = round(rp.FOV_corners(2) + ((0:1:(rp.searchGridPts-1))*scan_step(2)));
        % As our subimages are quite small, do preshifing of the source
        % images by the guess amounts, to reduce the amount that is thrown
        % away in alignment. Just using the same code as in aligned_FSC
        if ~isempty(param.guessx)
            switch sign(param.guessx)
                case 1
                    im0 = im0(:,1+param.guessx:end);
                    im1 = im1(:,1:end-param.guessx);
                case -1
                    im0 = im0(:,1:end+param.guessx);
                    im1 = im1(:,1-param.guessx:end);
            end
        end
        if ~isempty(param.guessy)
            switch sign(param.guessy)
                case 1
                    im0 = im0(1+param.guessy:end,:);
                    im1 = im1(1:end-param.guessy,:);
                case -1
                    im0 = im0(1:end+param.guessy,:);
                    im1 = im1(1-param.guessy:end,:);
            end
        end
        %preshifting is done, so set the guesses to 0 now.
        param.guessx = 0;
        param.guessy = 0;


        for ii = 1:rp.searchGridPts

            for jj = 1:rp.searchGridPts
                param.crop = {Rstarts(ii)+ (0:1:(rp.roisize(1) - 1)), ...
                              Cstarts(jj)+ (0:1:(rp.roisize(2) - 1))};
                rrr(ii,jj).mask(param.crop{:}) = 1;
                rrr(ii,jj).mask_centre = [(param.crop{1}(1) + param.crop{1}(end))/2,...
                                          (param.crop{2}(1) + param.crop{2}(end))/2];
                    
              [rrr(ii,jj).resolution, rrr(ii,jj).stat, rrr(ii,jj).FSC, rrr(ii,jj).T, rrr(ii,jj).freq, rrr(ii,jj).n] = ...
                  aligned_FSC(im0, im1, param);
            end
        end        
    end

    if rp.savefig
        ff = figure(100);  % aligned_FSC puts the main graph resolution graph in figure 100.               
        saveas(ff,rp.savefigname);
    end

    % if ~rp.autoplotclose
    %     close_figs = input('Enter c to close all figures, or any other key to keep them open.\n','s');
    %     if strcmpi(close_figs,'c')
    %         close('all');
    %     end
    % else

    if rp.autoplotclose
        pause(1);
        close('all');
    end

    if ~isempty(rp.saveresname)
        save(rp.saveresname, 'rrr','rp','param');
    end


end

function [img, im_dat] = load_check_rotate(fname)
    if ischar(fname)
        r0 = load(fname);
    elseif isstruct(fname)
        r0 = fname;
    end
    assert(isfield(r0,'recon_dat'),'Input contains no recon_dat field');
    if length(size(r0.recon_dat.obj_int{end})) == 2
        img = r0.recon_dat.obj_int{end};
    else
        img = r0.recon_dat.propogated_obj;
    end
    img = imrotate(img, -r0.recon_dat.rotation_angle);
    im_dat = r0.recon_dat;
end

function masked_im = make_and_mult_mask(im)
% in case we want to do masking at this stage:
masked_im = im;

end






