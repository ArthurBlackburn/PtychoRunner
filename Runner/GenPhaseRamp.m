% Produces a phase ramp required to produce a specified shift in the FT ('real-space') of
% the produced (complex) phase ramp.
%
% Usage:
%
% Ramp = GenPhaseRamp(RealShift, OutputSize)
%
%  RealShift - The number of pixels in real-space (i.e. the FT) the ramp produces [x, y].
%  OutputSize - The size of the required output matrix. [Columns 'x', Rows 'y']
%
%  Ramp - The output complex-numbered phase ramp matrix.
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

function Ramp = GenPhaseRamp(RealShift, OutputSize)
    
sub_x = ifftshift(-fix(OutputSize(1)/2):ceil(OutputSize(1)/2)-1)/OutputSize(1);
sub_y = ifftshift(-fix(OutputSize(2)/2):ceil(OutputSize(2)/2)-1)/OutputSize(2);

[sub_x,sub_y] = meshgrid(sub_x, sub_y);

Ramp = exp(1i*2*pi*(sub_x*RealShift(1) + sub_y*RealShift(2)));

