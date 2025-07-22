function [big_summary, obj_crp ] = savesummaryimage3( save_summary_image_name, obj, prb, ...
    original_rot_angle, description_ascii, options)
% Prepares a 'summary image' of a ptychographic phase image, preparing
% object and probe plots, and creating a FFT of complex image, etc.
% Top Left: Probe amplitude
% Top Middle: Object Phase (not unwrapped)
% Top Right: Object Amplitude
% Bottom Left: Probe Phase (not unwrapped)
% Bottom Middle: Ampltiude of the FT of the phase of the reconstruction.
% Bottom Right: Ampltiude of the FT of the amplitude of the reconstruction.
%
% Note: by default amplitude images have the extremal 1% of value clipped or saturated. This can
% cause some images to look unexpectedly saturated. For publications purposes you'll probably want
% to lower this value, which is held in the options structure in def_options.percent_sat to < 1, or
% use some image prep method.
%
% [big_summary, obj_crp ] = ...
%       savesummaryimage3( save_summary_image_name, obj, prb, ...
%                         original_rot_angle, description_ascii, options)
%
% Examples:
%
% savesummaryimage3(path_and_filename, ...
%   recon_dat.obj_int{end}, recon_dat.probe_int{end}, ...
%  -recon_dat.rotation_angle,...
%  'demo',struct('phase_pres_method','product'));
%
% [big_summary, obj_crp ] = savesummaryimage3('temp.tif', ...
%   recon_dat.propogated_obj, recon_dat.probe_int{end}, ...
%  -recon_dat.rotation_angle);
%
% savesummaryimage3(path_and_filename, ...
%   recon_dat.obj_int{end}, recon_dat.probe_int{end}, ...
%  -recon_dat.rotation_angle,...
%  'Some text that describes the reconstruction',...
%   struct('phase_pres_method','propa',...
%          'recon_pix_size', recon_dat.recon_pixel_size, ...
%          'lambda',recon_dat.lambda ,...
%          'Dz', recon_dat.recipe.params.delta_z_angstrom*1e-10));
% 
% options structure fields: {(default), other_vals}
%
%         give_windowed_fft: {(true),  false}
%         window_type: {('hann'), 'hamming'} 
%         phase_pres_method: {('sum_phases'), product}
%         save_sep: {(false),  true}
%         others are in the code... have a look!
%
% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2020
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
% 
% *******************************************************************

def_options.phase_pres_method = 'product';
def_options.save_sep = false;
def_options.window_type = 'hann';
def_options.give_windowed_fft = true;
def_options.edge_trim_DPSIZE = [];
def_options.recon_pix_size = [];
def_options.gamma = 0.4;
def_options.percent_sat = 1;

if nargin < 6
    options = def_options;
end

for fn = fieldnames(def_options)'
    if ~isfield(options, fn{1})
        options.(fn{1}) = def_options.(fn{1});
    end
end

phase_pres_method = options.phase_pres_method;
save_sep = options.save_sep;
window_type = options.window_type;
give_windowed_fft = options.give_windowed_fft;

if nargin < 5
    description_ascii = [];
end

obj_size = size(obj);
prb_size = size(prb);
if length(obj_size) > 2 && obj_size(end) > 1
    % For multislice summed phases of individual slices are likely to better
    % represent the object:

    switch phase_pres_method
        case 'sum_phase'
            obj_inds = num2cell(ones(1,length(obj_size)));
            obj_inds(1:2) = {':'};
            obj_phase = angle(squeeze(obj(obj_inds{:})));
            for slc = 2:obj_size(end)
                obj_inds{end} = slc;
                obj_phase = obj_phase + angle(squeeze(obj(obj_inds{:})));
            end
            obj = squeeze(prod(obj,length(obj_size)));
            % obj should be used for ampltiude which is still found using
            % product method. Phase and amplitude
            % are of course also correct in obj, but phase may more or less 
            % wraps in it compared to the summed phase of the individual slices...
        case 'product'
            obj_phase = angle( squeeze(prod(obj,length(obj_size))) );
            obj = squeeze(prod(obj,length(obj_size)));
        case 'propa' % short for propogate, and pun partly intended as this is a bit more proper.
            
            obj_info.recon_pixel_size = options.recon_pix_size;
            obj_info.lambda = options.lambda;
            obj_info.delta_z_m = options.Dz;
            if ~isfield(options,'wave_gen_func')
                options.wave_gen_func = [];                            
            end
            obj = propogate_object(obj, obj_info, options.wave_gen_func);
            obj_phase = angle(obj);

        otherwise        
            warning('Unknown phase presentation method, using default product instead');
            obj_phase = angle( squeeze(prod(obj,length(obj_size))) );
            obj = squeeze(prod(obj,length(obj_size)));
    end

else
    obj = squeeze(obj);
    obj_phase = angle(obj);
end


if length(prb_size) > 2
    if length(prb_size) == 3
        prb = squeeze(prb(:,:,1));
    elseif length(prb_size) == 4
        prb = squeeze(prb(:,:,1,1));
    else
        fprintf('Not really sure when dimensions of the probe are greater than 4, so summming along last dimension.\n');
        prb = squeeze(prod(prb,length(prb)));         
    end
end



crunch16  = @(x) uint16(double(intmax('uint16'))*mat2gray(x));
pres_data = @(x) imadjust(crunch16(x)); % use default percent sat of 1%
% 1% of pixels will be saturated by default.
pres_data_sat = @(x, percent_sat) imadjust (crunch16(x),...
    stretchlim(crunch16(x)), [percent_sat/100, 1-percent_sat/100]); % pre of pixels will be saturated by default.

if isempty(options.edge_trim_DPSIZE)
    DP_SIZE = 256; % don't know why I at one point hard-coded this....
else
    DP_SIZE = options.edge_trim_DPSIZE;
end
% DP_SIZE = size(prb,1);
over_factor = 1.5;
q = original_rot_angle;
obj_size = size(obj) - 2*ceil(over_factor*DP_SIZE/2);
h = obj_size(1);
w = obj_size(2);
rot_mat = [[sind(q) cosd(q)];[cosd(q) sind(q)]];
rot_mat = abs(rot_mat);
% inv_rot = inv(rot_mat);
% original_size = inv_rot*[w;h];
original_size = rot_mat\[w;h]; % same as above in Matlab lingo.
a = original_size(1);
b = original_size(2);

% Rotate the source image back, working with the abs and phs:
% rot_abs = imrotate(abs(obj), -q);
% rot_phs = imrotate(angle(obj), -q);

rot_abs = imrotate(abs(obj), q);
rot_phs = imrotate(obj_phase, q);


% Recombine and crop:
% rot_obj = rot_abs.*exp(2i*pi*rot_phs);

centre_point = floor(size(rot_abs)/2);
% centre_point = floor(size(rot_obj)/2);
use_rows = centre_point(1) - floor(b/2) + (0:1:(b-1));
use_cols = centre_point(2) - floor(a/2) + (0:1:(a-1));
obj_abs_crp = rot_abs(use_rows, use_cols);
obj_phs_crp = rot_phs(use_rows, use_cols);

% obj_crp = rot_obj(use_rows, use_cols);
obj_crp.abs = obj_abs_crp;
obj_crp.phs = obj_phs_crp;

% summary_image = obj_crp;
% summary_plot(:,:,1,1) = imadjust(uint16(double(intmax('uint16'))*mat2gray(angle(obj_crp))));
% summary_plot(:,:,1,2) = imadjust(uint16(double(intmax('uint16'))*mat2gray(abs(obj_crp))));
im_phs = pres_data_sat(obj_phs_crp,0.0);

% im_phs = pres_data_sat(angle(obj_crp),0.0);
% im_abs = pres_data(abs(obj_crp));
im_abs = pres_data_sat(obj_abs_crp, options.percent_sat);



my_im_montage = [im_phs, im_abs];

gamma_c = options.gamma;

% Y = fft2(obj_phs_crp);
% Nov 29, 2022:
Y = fft2(ifftshift(obj_phs_crp));
phsFT = abs(fftshift(Y)).^gamma_c;
fft_phs = pres_data_sat(phsFT, options.percent_sat);

% Y = fft2(obj_abs_crp.^2);
% Nov 29, 2022:
Y = fft2(ifftshift(obj_abs_crp.^2));
absFT = abs(fftshift(Y)).^gamma_c; 
fft_abs = pres_data_sat(absFT, options.percent_sat);
my_fft_montage = [fft_phs, fft_abs];

% and the probe images:
% get probe image size
prb_size = size(prb,1);
% sort out the scale factor:
display_scale_factor = size(im_phs,1)/prb_size;

if display_scale_factor > 1
    display_scale_factor = floor(display_scale_factor);
end

prb_im_abs = pres_data_sat(imresize(abs(prb),display_scale_factor), options.percent_sat);
prb_im_phs = crunch16(imresize(angle(prb),display_scale_factor));
black_pad = zeros([size(im_phs,1) - size(prb_im_phs,1),size(prb_im_phs,2)],'uint16');
prb_disp = [black_pad;...
            prb_im_abs;...
            prb_im_phs;...
            black_pad];


big_summary = [prb_disp, [my_im_montage; my_fft_montage]];

sz_ft_image = size(obj_crp.abs);
if ~isempty(options.recon_pix_size)
    one_over_d_pix_size_x_A = (1/(options.recon_pix_size*sz_ft_image(2)))/10^10;
    one_over_d_pix_size_y_A = (1/(options.recon_pix_size*sz_ft_image(1)))/10^10;
    one_over_pixsize = [one_over_d_pix_size_x_A, one_over_d_pix_size_y_A];
    generalscale_info_string = sprintf('Pixel Size: %3.2e, FFT pixel x: %3.2e, FFT pixel y:%3.2e',...
        options.recon_pix_size*10^10, one_over_pixsize);
    simple_desc = [description_ascii,', ',generalscale_info_string];
else
    general_info_string = [];
    simple_desc = description_ascii;
end

ArthurTiff(big_summary, save_summary_image_name, simple_desc);

if save_sep
    % output the individual abs and phase images without any suturation 
    % treatment etc, and squared abs of FFT of complex image.
    % The latter is simmilar to what would be seen in FFT
    % of the image in a TEM.
    [filepath,name,ext] = fileparts(save_summary_image_name);
    abs_im_name = fullfile(filepath,[name,'_abs',ext]);
    ArthurTiff(obj_crp.abs, abs_im_name, simple_desc);
    %
    phs_im_name = fullfile(filepath,[name,'_phs',ext]);
    ArthurTiff(obj_crp.phs, phs_im_name, simple_desc);
    %
    fft_desc = sprintf('%s, %s, gamma: %2.3f',description_ascii, generalscale_info_string, gamma_c);
    sname = fullfile(filepath,[name,'_absoffft_of_phs',ext]);
    ArthurTiff(fft_phs, sname, fft_desc);
    sname = fullfile(filepath,[name,'_absoffft_of_abssqd',ext]);
    ArthurTiff(fft_abs, sname, fft_desc);
    % also make uniformly scaled FFT, so round features should be round
    % etc:
    [fft_phs_equal, use_pixel_size, rescaled_image] = square_up(fft_phs);
    if rescaled_image
        equal_pix_size = sum(one_over_pixsize.*use_pixel_size);
        equal_info_string = sprintf('Pixel Size: %3.2e, FFT pixel: %3.2e',...
        options.recon_pix_size*10^10, equal_pix_size);
        equa_fft_desc = sprintf('%s, %s, gamma: %2.3e',description_ascii, equal_info_string, gamma_c);
        sname = fullfile(filepath,[name,'_absoffft_of_phs_equal',ext]);
        ArthurTiff(fft_phs_equal, sname, equa_fft_desc);
        
        fft_abs_equal = square_up(fft_abs);
        sname = fullfile(filepath,[name,'_absoffft_of_abssqd_equal',ext]);
        ArthurTiff(fft_abs_equal, sname, equa_fft_desc);
        
    end
    
    
    complex_obj_crp = obj_crp.abs .* exp(1i*obj_crp.phs);
    obj_crp.complex_obj_crp = complex_obj_crp;
    
    ft = fftshift(fft2(complex_obj_crp));
    ft_intens = abs(ft).^(2*gamma_c);
    fft_intens_name = fullfile(filepath,[name,'_abssqd_of_complexfft',ext]);
    ArthurTiff(ft_intens , fft_intens_name, simple_desc);
    % Maybe we also want to apply a hanning / hamming window before fft?
    if give_windowed_fft
        n_max = max(size(complex_obj_crp));
        if strcmpi(window_type,'hann')        
            window = .5*(1 - cos(2*pi*(0:n_max-1)'/(n_max-1)));
        end
        if strcmpi(window_type,'hamming')
            window = 0.54 - 0.46*cos(2*pi*(0:n_max-1)'/(n_max-1));
        end
        window = (window*window');
        window = imresize(window,size(complex_obj_crp));
        ftw = fftshift(fft2(window.*complex_obj_crp));
        ftw_intens = abs(ftw).^(2*gamma_c);
        ftw_intens_name = fullfile(filepath,[name,'_ftwintens',ext]);
        ArthurTiff(ftw_intens, ftw_intens_name, simple_desc);
    end
    
end

% for examples of color plotting see:
% C:\Users\Arthur\OneDrive\Documents\Papers\JSM2022\PaperFigures.m
% and 
% L:\em_utils\Basic_Utils\imagesc_hsv_ab.m

end

function [ft_image, use_pixel_size, rescaled_image] = square_up(ft_image)
    sz_ft_image = size(ft_image);
    if sz_ft_image(2) > sz_ft_image(1) % i.e. x_pixel in Fourier space is smaller
        % one_over_d_pix_size_A = one_over_d_pix_size_x_A;
        use_pixel_size = [1 0];
        ft_image = imresize(ft_image, [sz_ft_image(2), sz_ft_image(2)]); % make it square by streching in y;

        rescaled_image = true;
    elseif sz_ft_image(2) < sz_ft_image(1)

        % one_over_d_pix_size_A = one_over_d_pix_size_y_A;
        use_pixel_size = [0 1];
        ft_image = imresize(ft_image, [sz_ft_image(1), sz_ft_image(1)]); % make it square by streching in y;
        rescaled_image = true;

    else
        rescaled_image = false;
        use_pixel_size = [1 0];
        % one_over_d_pix_size_A = one_over_d_pix_size_x_A;
    end
end
