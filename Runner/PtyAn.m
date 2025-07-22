classdef PtyAn
% PtyAn: Class of convenience, bundling together methods for ptychography analysis. 
% Contains things such as SNR analysis, PCTF extraction etc. Eventually this may grow 
% into a something bigger and better for analyzing a ptychographic reconstructions.
%
% For examples of usage see "...\TestScripts\Part_5_PCTF_SNR_on_SyntheticData.m"
%
% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2024
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
% 
% *******************************************************************
    
    properties
        Property1
    end
    
    methods(Static)
        
        function gt = get_trimmed(gt, active_pixels)
            gt.centre = ceil(size(gt.object)/2);
            gt.start_ind = gt.centre - floor(active_pixels/2);
            gt.central_trimmed = gt.object(gt.start_ind(1)+(1:active_pixels), gt.start_ind(2)+(1:active_pixels));
        end

        function plot_quick(gt, rc, field_name)
            figure;
            a1 = subplot(2,2,1);
            imagesc(abs(gt.(field_name))); axis image; title('Abs GT');
            a2 = subplot(2,2,2);
            imagesc(angle(gt.(field_name))); axis image; title('Angle GT');
            a3 = subplot(2,2,3);
            imagesc(abs(rc.(field_name))); axis image; title('Abs RC');
            a4 = subplot(2,2,4);
            imagesc(angle(rc.(field_name))); axis image; title('Angle RC');
            sgtitle(field_name,'Interpreter','none');

            linkaxes([a1 a2 a3 a4]);
        end

        function rcn = get_TFs(gt, rcn, real_space_pix_size, kwargs)

            % gt. is the ground truth, rc.(rc_name) is reconstructed object. Both are complex
            % images. Probably best to use rc_name = 'best_match' here, as this has
            % not already been normed etc (unlike best_match_normed).

            if nargin < 3
                real_space_pix_size = 1;
            end
            % PCTF
            % 1) Simple approach

            p = inputParser;
            p.addParameter('rc_name','best_match',@(x) ismember(x,{'best_match','best_match_normed','best_match_levelled'}));
            p.addParameter('tlevel',0.35, @(x) isnumeric(x));  % for threshholded pctf this is the threshold..
            p.addParameter('tftype','pctf',@(x) ismember(x,{'pctf','actf'}));
            p.addParameter('nbins',50, @(x) isnumeric(x));  % number of radial bins up to extent of fft
            p.addParameter('for_correction',false, @(x) islogical(x));  % whether to create a 2d map
            p.parse(kwargs{:});
            s = p.Results;
            
            
            switch s.tftype
                case 'pctf'
                    % get the phase unwrapped first:
                    obj_meas = Unwrap_TIE_DCT_Iter(angle(rcn.(s.rc_name)),60*4);
                    gt_meas = Unwrap_TIE_DCT_Iter(angle(gt),60*4);
                otherwise
                    error('Invalid ctf type\n');
            end

            obj_meas = obj_meas - mean(obj_meas(:));
            gt_meas = gt_meas - mean(gt_meas(:));
            
            pc_simple = abs(ft2(obj_meas,1)) ./ abs(ft2(gt_meas,1));
            rcn.([s.tftype,'_simple']) = PtyAn.get_radial_aves_simple(pc_simple, s.nbins, true(size(pc_simple)), s.for_correction);
            rcn.([s.tftype,'_simple']).radial_ord =  rcn.([s.tftype,'_simple']).radial_ord * (1/(real_space_pix_size*size(pc_simple,1)));
            rcn.([s.tftype,'_4bin_simple']) = PtyAn.get_radial_aves_simple(PtyAn.bin_image(pc_simple), s.nbins, [], s.for_correction);          
            % 2) Thresholded approach: ignore phase measurements coming from
            % frequency components that have relatively small intensity in the
            % image. First look at the intensities (ampl.^2), find which regions
            % are strong, then compare these.
            abs_fft_ampsqd = abs(ft2(abs(rcn.(s.rc_name)).^2,1));
            filtered_fft = stdfilt(abs_fft_ampsqd,true(5)).^0.4;
            % somewhat arbitrary but going for a simple thresh on this:
            thresh_level = s.tlevel * max(filtered_fft,[],'all');
            fft_mask = filtered_fft > thresh_level;
            rcn.([(s.tftype),'_threshed']) = PtyAn.get_radial_aves_simple(pc_simple, s.nbins, fft_mask, s.for_correction);
            rcn.([s.tftype,'_threshed']).radial_ord =  rcn.([s.tftype,'_threshed']).radial_ord * (1/(real_space_pix_size*size(pc_simple,1)));


            rcn.(s.tftype).simple_2d = pc_simple;            
            rcn.(s.tftype).fft_mask = fft_mask;
            rcn.(s.tftype).thresh_level = thresh_level;
            rcn.(s.tftype).abs_fft_ampsqd = abs_fft_ampsqd;
            rcn.(s.tftype).unwrapped_obj = obj_meas;
            rcn.(s.tftype).unwrapped_gt = gt_meas;

        end
        

        function fhandle = plot_TFs(rc, poss_fields, fhandle, ymodifier, leg_cell)
            if nargin < 2
                poss_fields = {'pctf_simple','pctf_threshed'};    
            end
            if nargin < 3
                fhandle = [];
            end
            if nargin < 4
                ymodifier = @(y) y;
            end
            if nargin < 5
                leg_cell = {};
            end
            available_fields = fieldnames(rc);
            fmask = false(length(available_fields),1);
            for nn = 1:length(poss_fields)
                fmask = fmask | strcmpi(available_fields, poss_fields{nn});
            end
            num_plots = sum(fmask,'all');
            if num_plots == 0
                fprintf('plot_TFs() has nothing to plot\n');
                return; % i.e. don't go any further.
            end
            plot_fields = available_fields(fmask);
            if isempty(leg_cell)
                leg_cell = plot_fields;
            else

            end

            if isempty(fhandle)
                fhandle = figure;
            else
                figure(fhandle)
            end
            
            for plt_ind = 1:num_plots
                plot(rc.(plot_fields{plt_ind}).radial_ord(2:end), ymodifier(rc.(plot_fields{plt_ind}).radial_ave(2:end)),'-x');
                hold on;    
            end
            legend(leg_cell,'Interpreter','none');

        end


        function aves = get_radial_aves_simple(im, nbnds, extra_mask, create_2d_map, dont_use_first)
            
            if nargin < 5
                % Sometimes zeroth point of the radial avarages we are calculating here can be Inf, 
                % if zero-frequency components are set equal (giving noise = zero), so you
                % might want to exclude it.
                dont_use_first = 0;
            end
                                       
            if nargin < 4
                create_2d_map = false;
            end
                
            if nargin < 3 || isempty(extra_mask)
                % extra_mask defined such that if there are any data point you
                % want excluded from the average, set these in the mask to
                % be zero. 
                extra_mask = true(size(im));
            end
    
            centre_index = @(nr)(nr-rem(nr,2))/2+1;
            centre_pt = centre_index(size(im));

            [nrows, ncols] = size(im);
    
            [xx, yy] = meshgrid((1:ncols) - centre_pt(2), (1:nrows) - centre_pt(1));
            rr = sqrt(xx.^2 + yy.^2);
            maxr = min([centre_pt - 1, size(im) - centre_pt]);
            bnd_step = maxr/nbnds;
            aves.radial_ave = zeros(1,nbnds);
            if create_2d_map
                aves.map2d = zeros(size(im));
            end
    
            for ii = 1:nbnds
                mask = ((ii - 1)*bnd_step <= rr) & (rr < (ii*bnd_step));
                mask = extra_mask & mask;
                aves.radial_ave(ii) = sum(im(mask),'all') / sum(mask,'all');
            end
            
            if create_2d_map
                for ii = 1:nbnds
                    mask = ((ii - 1)*bnd_step <= rr) & (rr < (ii*bnd_step));
                    mask = extra_mask & mask;
                    if ii > dont_use_first
                        aves.map2d(mask) = aves.radial_ave(ii);
                    else
                        aves.map2d(mask) = aves.radial_ave( dont_use_first + 1 );
                    end
                end
            end            
    
            aves.radial_ord = bnd_step*(1:1:nbnds);
            % remove the elements that are 0/0:
            aves.radial_ord = aves.radial_ord(~isnan(aves.radial_ave));
            aves.radial_ave = aves.radial_ave(~isnan(aves.radial_ave));

        end

        function rc = get_SNR_measures(rc, gt, real_space_pix_size, kwargs)
            
            if nargin < 3
                real_space_pix_size = 1;
            end

            p = inputParser;            
            p.addParameter('amp_mult',1, @(x) isnumeric(x)); % if ampltiude needs scaling... 
            %  probe mode and probe contraint issues, like removing edge
            %  counts lead to some intensity changes, but overall it is
            %  usually the changes in amplitude that matter as the signal,
            %  so amplitude scaling seems acceptable.
            p.addParameter('phase_mult',1, @(x) isnumeric(x)); % if you want to see multipliers of phase do,
            % then here's your chance. Can't think it would be necessary.
            % Physics of phase is robust!
            p.addParameter('phase_shift',0, @(x) isnumeric(x)); % phase shift might be needed due 
            % to arbitray assignment of mean phase in both object and reconstruction.
            p.addParameter('gt_fieldname','best_match', @(x) ischar(x));
            p.addParameter('nbins',50, @(x) isnumeric(x));  % number of radial bins up to extent of fft
            p.addParameter('for_correction',false, @(x) islogical(x));  % whether to create a 2d map
            p.addParameter('recon_fieldname','best_match', @(x) ischar(x));
            p.addParameter('compon', 'complex'); % if you want to look at a component of the signal / noise, specify it 
                                              % here using a function eg  @(x) angle(x)
            p.addParameter('bin_factor',4);
            p.parse(kwargs{:});
            s = p.Results;
            
            if strcmp(s.compon,'complex')
                   % The 'signal' (which is also our ground truth)
                   signal = gt.(s.gt_fieldname);
                   % Determine the noise
                   % a) the object:
                   obj = rc.(s.recon_fieldname);                   
                   % b) adjust the output as appropriate:                   
                   modified_rc = s.amp_mult*abs(obj) .* exp(1i*(angle(obj)*s.phase_mult + s.phase_shift) );
                   % finally get the (complex) noise 

                   noise = modified_rc - signal;
            elseif strcmpi(s.compon, 'angle')
                   signal = angle(gt.(s.gt_fieldname));
                   obj = angle(rc.(s.recon_fieldname));
                   modified_rc = s.phase_mult*obj + s.phase_shift;
                   noise = modified_rc - signal;
            elseif strcmpi(s.compon, 'angle_flatmod')
                   signal = exp(1i*angle(gt.(s.gt_fieldname)));
                   obj = exp(1i*angle(rc.(s.recon_fieldname)));
                   modified_rc = exp(1i*s.phase_mult*angle(obj) + s.phase_shift);
                   noise = modified_rc - signal;
            else
                error('unrecognized measurement option.')
            end

            % Get the FTs:
            % NOTE In tests performed here, we _purposely choose a square
            % image_ so thet size in x is same in y, so x and y frequency
            % range is the same by fft. If that's not the case, we'd need
            % to do some interpolation and scaling to make sure we're
            % really getting radial average in the above. The FRC / FSC
            % implements this more advanced approach. However, here in
            % general we have the spatial frequency per pixel in the FTs
            % as:
            freq_per_pixel = (s.bin_factor/(real_space_pix_size*size(modified_rc,1)));

            rc.spatial_freq.row_1axis = (1/(real_space_pix_size*size(modified_rc,1)))*...
                (-floor(size(modified_rc,1)/2):1:ceil(size(modified_rc,1)/2)-1);
            rc.spatial_freq.col_2axis = (1/(real_space_pix_size*size(modified_rc,2)))*...
                (-floor(size(modified_rc,2)/2):1:ceil(size(modified_rc,2)/2)-1); 

            % Get signal and noise and their radial averages from
            % Thibault P and Guizar-Sicairos M (2012) Maximum-likelihood refinement for coherent 
            % diffractive imaging. New Journal of Physics 14, 063004.
            % Equation 43
            rc.S = PtyAn.bin_image(abs(ft2( signal ,1).^2), s.bin_factor);

            rc.N = PtyAn.bin_image(abs(ft2( noise  ,1)).^2, s.bin_factor);
            
            rc.S_rad = PtyAn.get_radial_aves_simple(abs(rc.S), s.nbins, [], s.for_correction);
            rc.S_rad.log10_radial_ave = log10(rc.S_rad.radial_ave);

            rc.N_rad = PtyAn.get_radial_aves_simple(abs(rc.N), s.nbins, [], s.for_correction);
            rc.N_rad.log10_radial_ave = log10(rc.N_rad.radial_ave);

            % If there is no noise (ie. N =0) at certain points in the FFT, we 
            % will want to be aware of this. It also raises a bit of
            % dilemma on how to calculate SNR average, as if one point in
            % the series of SNR is Infinite ( arising from N = 0), then the
            % average of all the points is inifinite. This would be
            % be a misleading measure, imo, as all the other points would have 
            % reasonable value. In experimental data, the probability of
            % this probably effectivesly zero, but in 'fake' infinite dose
            % data, it's a maybe...
            if any(rc.N == 0,'all') 
                warning('Calculated noise contains zeros. Carefully check result validity.');
                % There might be (almost is sure to be) better ways of
                % dealing with this case, but for now I'll try:
                rc.N(n_zero) = min(N(~n_zero),[],'all'); 
                % if the warning above does not appear, you don't have
                % worry about whether the above is fair or not... but the
                % above just replaces zeros values the next biggest non-zero number
                % in the set.
            end

            rc.SNR = rc.S./rc.N;
            
            rc.SNR_rad = PtyAn.get_radial_aves_simple(abs(rc.SNR), s.nbins, [], s.for_correction);
            rc.SNR_rad.log10_radial_ave = log10(rc.SNR_rad.radial_ave);
            rc.SNR_rad.radial_ord = freq_per_pixel * rc.SNR_rad.radial_ord;

            rc.SNR_rad_lin = PtyAn.get_radial_aves_simple(abs(rc.SNR), s.nbins, [], s.for_correction);

            rc.phase_mult = s.phase_mult;
            rc.amp_mult = s.amp_mult;
            rc.phase_shift = s.phase_shift;

        end

        function im = bin_image(im, bin_factor)
            if nargin < 2
                bin_factor = 4;
            end
            im = imresize(im, 1/bin_factor, 'bilinear');
        end

        function aa = plot_SNR(rc, aa, subtype)
            if isempty(aa)
                figure;
                a1 = subplot(1,3,1);
                a2 = subplot(1,3,2);
                a3 = subplot(1,3,3);
                aa = [a1, a2, a3];
            end
            if nargin < 3
                subtype = [];
            end
            axes(aa(1));
            plot(rc.(['SNR',subtype,'_rad']).radial_ave); hold on;
            title(['SNR',subtype]);
            axes(aa(2));
            plot(rc.(['S',subtype,'_rad']).radial_ave); hold on;
            title(['S',subtype]);
            axes(aa(3));
            plot(rc.(['N',subtype,'_rad']).radial_ave); hold on;
            title(['N',subtype]);
        end

        
        function [rc, fit_res] = phase_level(rc, reference_im, src_im, levelled_fieldname)

            % Adjusts the linear phase ramp background on a reconstruction,
            % so that it matches the linear phase ramp background of the
            % ground truth. Small phase ramps accross the entire image have
            % little / no impact on the overall useful information in the
            % image and can be removed. Sometimes appear in reconstruction,
            % if the solver stumbles into tilting the entire object phase
            % rather than just the probe phase, in order to get the forward
            % propogated wave FFT (particularly centre of the bright field
            % disc) to match the experimental data. Might be an idea for
            % future work to periodically check for phase ramp in object,
            % and periodically push that global tilt on to the probe... Or
            % did someone do that already within 'tilt corrected
            % ptycho'?.. Hmmm.

            ref_unwrap = Unwrap_TIE_DCT_Iter(angle(reference_im),60*4);
            src_unwrap = Unwrap_TIE_DCT_Iter(angle(src_im),60*4);
            [xx, yy] = meshgrid(((1:size(ref_unwrap,2))-1)/ (size(ref_unwrap,2)-1),...
                                ((1:size(ref_unwrap,1))-1)/ (size(ref_unwrap,1)-1));
            [xData, yData, zData_ref] = prepareSurfaceData( xx, yy, double(ref_unwrap));
            ft = fittype( 'poly11' );
            [fit_res.ref_model, fit_res.gof_ref] = fit( [xData, yData], zData_ref, ft );
            [xData, yData, zData_src] = prepareSurfaceData( xx, yy, double(src_unwrap));
            ft = fittype( 'poly11' );
            [fit_res.src_model,fit_res.gof_src] = fit( [xData, yData], zData_src, ft );
            % make a version of the src, that has same background as src:
            src_corrected = src_unwrap - fit_res.src_model(xx,yy) + fit_res.ref_model(xx,yy);
            rc.(levelled_fieldname) = abs(src_im) .* exp(1i*src_corrected);            
            
        end
        
        function [rc, gt, best_ind] = simple_align(rc, gt, bi, measure, kwargs)

            % Brute force alignment search. Covers the entire space
            % specified. With fake images with many feature a more
            % optimized approach could be taken, but with lattice images of
            % common in materials related HRTEM for example, there are many 
            % 1000's of local minima and other quick search routines thus need 
            % careful babysitting or filtering particularly in presence of
            % noise. A full search just works but is compuationally a
            % little expensive.

            if nargin < 3
                bi = [];
            end

            if nargin < 4
                 measure = @(x) angle(x);
            end
            
            p = inputParser;
            p.addParameter('width_spec','ratio',@(x) ismember(x,{'ratio','single'}));
            p.addParameter('offset_extent',30, @(x) isnumeric(x));
            p.addParameter('search_step',1, @(x) isnumeric(x))
            p.addParameter('subsample',1, @(x) isnumeric(x));
            p.addParameter('match_field','central_trimmed', @(x) ischar(x));
            p.parse(kwargs{:});
            s = p.Results;

            search_range = 1:s.search_step:(2*s.offset_extent + 1);

            if isempty(bi)
                xcms = zeros(numel(search_range), numel(search_range));
                aa = measure(gt.(s.match_field)((s.offset_extent+1):s.subsample:(end-s.offset_extent),...
                                       (s.offset_extent+1):s.subsample:(end-s.offset_extent)));

                for iii = 1:1:numel(search_range)
                    ii = search_range(iii);
                    for jjj = 1:1:numel(search_range)
                        jj = search_range(jjj);
                        bb = measure(rc.(s.match_field)(ii:s.subsample:(ii + end - search_range(end)),...
                                                        jj:s.subsample:(jj + end - search_range(end))));
                        xcms(iii, jjj) = sum(sum(aa.*bb));
                    end
                end

                figure;
                imagesc(xcms); colorbar; axis image;
                title('Cross Correlation Measure vs. X & Y Offset')

                [~, m_ind] = max(xcms,[],'all','linear');
                [r_max_ind,c_max_ind] = ind2sub(size(xcms),m_ind);
                r_max = search_range(r_max_ind);
                c_max = search_range(c_max_ind);
            else
                r_max = bi(1);
                c_max = bi(2);
            end
            rc.best_match = rc.(s.match_field)(r_max:(r_max + end - search_range(end)),...
                                               c_max:(c_max + end - search_range(end)));
            gt.best_match = gt.(s.match_field)((s.offset_extent+1):(end-s.offset_extent),...
                                               (s.offset_extent+1):(end-s.offset_extent));
            best_ind = [r_max,c_max];
        end        
        
    end
end

