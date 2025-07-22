% Moves the specified point of an image to the centre using Fourier shift theorem. 
%
% Usage:
% im_out = Recentre(im_in, new_centre_xy) 
% 
% Inputs
%  im_in - The image to be shifted.
%  nnew_centre_xy - Co-ordinates to move to the centre, 
%                   [cols 'x', row 'y'] (can be fractional)
%
% Outputs
%  im_out - The recentred image image.
%
% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2018
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
% 
% *******************************************************************

function im_out = Recentre(im_in, new_centre_xy)
    % The current centre of the input image.
    % ipCentre = round(size(im_in)/2) + 1; % this is one possibility.
    % But if we're dealing with a FT, it is more appropriate to use:
    ipCentre = floor(size(im_in)/2) + 1;    % [row 'y', col 'x'];
    ipCentre_xy = ipCentre(end:-1:1);       % [col 'x', row 'y'];
        
    % Recentre using Fourier Shift theorem
    displacement_xy = new_centre_xy - ipCentre_xy;
    im_size_xy = size(im_in);
    im_size_xy = im_size_xy(end:-1:1);

    ramp = GenPhaseRamp(displacement_xy, im_size_xy);
    ftInp = fft2(im_in);
    im_out = ifft2(ftInp.*ramp);
end
