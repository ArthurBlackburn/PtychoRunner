function G = ft2(g, delta)
% A simple function for finding the 2D FFT of function, with scaling for the step size (delta)
% usage:
% G = ft2(g, delta)
% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2021
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
% 
% *******************************************************************

G = fftshift(fft2(fftshift(g))) * delta^2;