function im1 = slice_progate(im0, pix_size, lambda, Dz)
% Takes a source image wave in and propogates it, padding it if necessary
% as a first step, before propogating then cropping the image back to same
% size as the output. Called it progate, because propagate was just a bit
% too proper (propa :-)).
%
% Usage:
%
% im1 = slice_progate(source_image, pix_size, lambda, delta_z)
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

pad_factor = 2;
im_shape = size(im0);
padded_sq_size = round(pad_factor * max(im_shape));
padded_sq_size = padded_sq_size + mod(padded_sq_size,2); % make it even

row_start = (padded_sq_size/2) - round(im_shape(1)/2);
row_end = row_start + im_shape(1) - 1;
col_start = (padded_sq_size/2) - round(im_shape(2)/2);
col_end = col_start + im_shape(2) - 1;

padded_im0 = ones(padded_sq_size,'like',im0);
padded_im0(row_start:row_end, col_start:col_end) = im0;

padded_projected_wave = Fresnel_fwd_propogate(padded_im0, pix_size, lambda, Dz);

im1 = padded_projected_wave(row_start:row_end, col_start:col_end);

end

function projected_wave = Fresnel_fwd_propogate(source_wave, pix_size, lambda, Dz)
    % Function to determine the propagated wave, when the propagtion distance
    % is small, i.e. Fresnel approach should be used (not farfield)
    %
    %  projected_wave = Fresnel_fwd_propogate(source_wave, pix_size, lambda, Dz)
    %
    %  (Source wave is assumed to be square)
    
    [~, ~, projected_wave] = ang_spec_prop(source_wave, lambda, pix_size, pix_size, Dz);
end


function [x2, y2, Uout] = ang_spec_prop(Uin, wvl, d1, d2, Dz)
    % Angular spectrum method of evaluating Fresnel diffraction integral, as 
    % applied by J. D. Schimdt
    %
    % function [x2 y2 Uout] = ang_spec_prop(Uin, wvl, d1, d2, Dz)
    
    N = size(Uin,1);   % assume square grid
    k = 2*pi/wvl;      % optical wavevector
    % source-plane coordinates
    [x1, y1] = meshgrid((-N/2 : 1 : N/2 - 1) * d1);
    r1sq = x1.^2 + y1.^2;
    % spatial frequencies (of source plane)
    df1 = 1 / (N*d1);
    [fX, fY] = meshgrid((-N/2 : 1 : N/2 - 1) * df1);
    fsq = fX.^2 + fY.^2;
    % scaling parameter
    m = d2/d1;
    % observation-plane coordinates
    [x2, y2] = meshgrid((-N/2 : 1 : N/2 - 1) * d2);
    r2sq = x2.^2 + y2.^2;
    % quadratic phase factors
    Q1 = exp(1i*k/2*(1-m)/Dz*r1sq);
    Q2 = exp(-1i*pi^2*2*Dz/m/k*fsq);
    Q3 = exp(1i*k/2*(m-1)/(m*Dz)*r2sq);
    % compute the propagated field
    Uout = Q3.* ift2(Q2 .* ft2(Q1 .* Uin / m, d1), df1);

end


