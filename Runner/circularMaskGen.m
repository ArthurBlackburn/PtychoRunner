function maskout = circularMaskGen(DPSIZE, circular_mask_radius)
% A very simple function for creating a circular 'mask' for a particular Fourier transform size,
% DPSIZE. The mask is centred at the zero point of (assumed square DP).
% 
% Usage: 
% mask = circularMaskGen(DPSIZE, circular_mask_radius)
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
    xrng = floor(-DPSIZE/2):1:((DPSIZE/2)-1);
    [XX, YY] = meshgrid(xrng);
    RR = sqrt(XX.^2 + YY.^2);
    maskout = RR <= circular_mask_radius;

end