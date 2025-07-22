function [res, ft] = radial_FT_from_recon(recon_dat, varargin)
% Used to determine the radial average of a reconstruction, 
% Some examples:
% e.g.
% 
%     ft_prof = radial_FT_from_recon(r0.recon_dat, 'do_plots',false, 'ave_meth','mean_profs');
%     ft_hist = radial_FT_from_recon(r0.recon_dat, 'do_plots',false, 'ave_meth', 'sums');
%     figure; 
%     plot(ft_prof.radial_ord(50:end), ...
%          ft_prof.radial_ave(50:end)./(max(ft_prof.radial_ave(50:end),[],'all')),'r');
%     hold on; 
%     plot(ft_hist.radial_ord(15:end), ...
%          ft_hist.radial_ave(15:end)./(max(ft_hist.radial_ave(15:end),[],'all')),'b');
%     legend({'average profiles','radial sums'});
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


    p = inputParser;
    p.addParameter('window_type','hamming', @(x)ismember(lower(x), {'hamming', 'hann'}));
    % In case you want to trim more from the edges of the images before taking FT.
    p.addParameter('custom_edgesize',[]);
    p.addParameter('do_plots',false);
    p.addParameter('ave_meth','sums',@(x)ismember(lower(x), {'sums', 'mean_profs'}));
    p.addParameter('n_pts',500);
    if isstruct(recon_dat)
        p.addParameter('lambda',recon_dat.lambda);
        p.addParameter('recon_pix_size',recon_dat.recon_pixel_size);    
    elseif ismatrix(recon_dat)
        p.addParameter('lambda',1);
        p.addParameter('recon_pix_size',1);    
    end

    p.parse(varargin{:});

    params = p.Results;

    if isstruct(recon_dat)
        if isfield(recon_dat,'propogated_obj')
            ob_in = recon_dat.propogated_obj;
            pres_meth = 'product';
        else
            ob_in = recon_dat.obj_int{1};
            pres_meth = 'propa';
        end
        % use save summary3 to get the cropped image:
        [~, obj_crp ] = savesummaryimage3([], ...
          ob_in, recon_dat.probe_int{end}, ...
         -recon_dat.rotation_angle,'test',...
         struct('edge_trim_DPSIZE',params.custom_edgesize,...
         'phase_pres_method',pres_meth));
        complex_obj_crp = obj_crp.abs .* exp(1i*obj_crp.phs);
    else
        complex_obj_crp = recon_dat;
    end

    n_max = max(size(complex_obj_crp));
    if strcmpi(params.window_type,'hann')        
        window = .5*(1 - cos(2*pi*(0:n_max-1)'/(n_max-1)));
    end
    if strcmpi(params.window_type,'hamming')
        window = 0.54 - 0.46*cos(2*pi*(0:n_max-1)'/(n_max-1));
    end

    window = (window*window');
    window = imresize(window,size(complex_obj_crp));

    % The WINDOWED FT:
    ft.win = fftshift(fft2(fftshift(window.*complex_obj_crp)));
    % The NON-WINDOWED ('raw') FT:
    ft.raw = fftshift(fft2(fftshift(complex_obj_crp)));
    % The image that we present of the (complex) FT:
    ft.image = abs(ft.win).^2;
    

    sz_ft_src_image = size(ft.image);    
    if sz_ft_src_image(1) ~= sz_ft_src_image(2)
        warning('Source image is not square, so fft will be resized using bicubic interpolation');
    end

    centre_index = @(nr)((nr-rem(nr,2))/2)+1;
    centre_point_RC_src = centre_index(size(ft.image));

    ft.src.one_over_d_pix_size_x_A = (1/(params.recon_pix_size*sz_ft_src_image(2)))/10^10;
    ft.src.one_over_d_pix_size_y_A = (1/(params.recon_pix_size*sz_ft_src_image(1)))/10^10;
    

    if sz_ft_src_image(1) ~= sz_ft_src_image(2)
        % Rescale the FTs, but save the source images, in case wanted, before scaling:
        ft.src.win = ft.win; ft.src.raw = ft.raw; ft.src.image = ft.image;
        ft.src.centre_point = centre_point_RC_src;

        if sz_ft_src_image(2) > sz_ft_src_image(1) % i.e. x_pixel in Fourier space is smaller        
            ft.one_over_d_pix_size_A = ft.src.one_over_d_pix_size_x_A;
            use_index = 2;
        else
            ft.one_over_d_pix_size_A = ft.src.one_over_d_pix_size_y_A;
            use_index = 1;
        end

        rescaled_image = 1;

        outsize = [sz_ft_src_image(use_index), sz_ft_src_image(use_index)];
        ft.win =   imresize(ft.win,   outsize); 
        ft.raw =   imresize(ft.raw,   outsize); 
        ft.image = imresize(ft.image, outsize); % make it square by streching in y (rows);
        ft.centre_point = 1 + round (...
            outsize./sz_ft_src_image .* ...
            (centre_point_RC_src - 1) );
        
    else
        rescaled_image = 0;
        ft.one_over_d_pix_size_A = ft.src.one_over_d_pix_size_x_A;
        ft.centre_point = centre_point_RC_src;
    end


    if params.do_plots(1)
        figure;
        subplot(1,1+rescaled_image,1);
        imagesc(ft.src.image, math.sp_quantile(ft.src.image,[1e-2, 1-1e-2],10) );
        title(sprintf('Abs(FT of source image).^2: size [Rows (Y) %d, Cols (X) %d]', size(ft.src.image)),'Interpreter','none');
        axis image; 
        colorbar
        hold on;
        plot(ft.src.centre_point(2), ft.src.centre_point(1), 'rx', 'MarkerSize', 15);
    
    
        if rescaled_image
            subplot(1,2,2);
            imagesc(ft.image, math.sp_quantile(ft.image,[1e-2, 1-1e-2],10) );
            title(sprintf('Abs(Equi-Scaled FT of source image).^2: size [Rows (Y) %d, Cols (X) %d]', size(ft.image)),'Interpreter','none');
            axis image; 
            colorbar
            hold on;
            plot(ft.centre_point(2), ft.centre_point(1), 'rx', 'MarkerSize', 15);
        end
    end
    % Just an in-passing note: For experimental reconstructions, there is some chance of phase wedge
    % existing in the (complex reconstruction). For most experiments this shoudl be quite small, and
    % hence the direct beam (Centre of the diffraction pattern for radial integration purposes) is
    % the same as the FT centre. However, if the wedge is large this would not be the case. Kind of
    % an edge case, but might neeed to think about adding a check for this here, to see if it is
    % significant.

    ft.obj = complex_obj_crp;
    res = radial_prof( ft.image, ft.centre_point, params.ave_meth, params.n_pts, ft.one_over_d_pix_size_A);


end


function res = radial_prof(src_mat, centre_pt, meth, npts, pix_size)
    
    % Use histogram type method (which uses all pixels), or average interpolated radial profiles?
    % If the data is smooth & non-sparse there might be a speed up with mean profiles method?..
    % Another thing to investigate...
    if nargin < 3
        meth = 'sums';
    end

    n_rows = size(src_mat,1);
    n_cols = size(src_mat,2);
    % I am averaging only over the radial range in which we have complete 360 degree coverage here.
    % Could remove, but for now I like the idea of radial avearge coming for a complete population,
    % rather than an a incomplete sample. For ptycho that's fine, as most of the stuff near edges is
    % noise anyway.

    row_range = min(n_rows-centre_pt(1), centre_pt(1)-1);
    col_range = min(n_cols-centre_pt(2), centre_pt(2)-1);
    
    radial_range = min(row_range, col_range);

    cropped_mat = src_mat( (-radial_range:1:radial_range) + centre_pt(1), ...
                           (-radial_range:1:radial_range) + centre_pt(2));

    % the interpolation coords:
    [X, Y] = meshgrid(-radial_range:1:radial_range);


    tic;
    if strcmpi(meth,'sums')

        if isempty(npts)
            npts = radial_range; % the full range;
        end
        radial_bands = linspace(1,radial_range,npts);
        res.radial_ave = zeros(size(radial_bands));
 
        % xor = @(P,Q) (P|Q)&~(P&Q); xor is built in...

        R = sqrt(X.^2 + Y.^2);
        P = false(size(R));

        for nn = 1:npts
            Q = R <= radial_bands(nn);
            mask_cur = xor(P,Q);
            P = Q;
            res.radial_ave(nn) = sum(cropped_mat(mask_cur),'all');
        end

        % adjust pixel size:
        % put the points at the centres of the bins
        res.radial_ord = linspace(pix_size,radial_range*pix_size,npts) - (pix_size/2); 

    elseif strcmpi(meth,'mean_profs')

        number_of_intervals = 180;  % how many radial profiles to average over?
        theta_step = 360/number_of_intervals; 
        radial_span = 0:1:radial_range;        
        theta_range = (2*pi/360) * (0:theta_step:(360-theta_step));
    
        % create the polar grid:
        v_x = (radial_span'*cos(theta_range))';
        v_y = (radial_span'*sin(theta_range))';

        radial_pats = interp2(...
            X, Y, cropped_mat , ...
            v_x, v_y, 'spline');

        res.radial_ave = mean(radial_pats,1);
        res.radial_ord = radial_span * pix_size;

    else
        error('unknown method');
    end



    fprintf('Radial Average Computed in %3.2f sec\n',toc);

end