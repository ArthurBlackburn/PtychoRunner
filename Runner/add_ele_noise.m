function [image_out, stats] = add_ele_noise( image_in )
% Take an image in where the pixel values represent the mean number of electron counts in that pixel,
% and then adds to it some random Poisson noise (based on the electron counts). If you are looking 
% For repeatable results the Matlab random number generator would need to be re-seeded prior to 
% calling this function. 
%
% The function also produces some statistics on the resulting image passed out.
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


% 1e12 / 1e-12 is the multiplier used by Matlab for double type number. See Matlab documentation
image_out = round(1e12*imnoise(1e-12*double(image_in),'poisson'));

if nargout == 2
    stats.mean = mean(double(image_out(:)));
    stats.rootmean = sqrt(stats.mean);
    stats.stddev_out_minus_in = std( double(image_out(:)) - double(image_in(:)));
end

end

