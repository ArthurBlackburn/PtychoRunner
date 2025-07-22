function [circshifted_dat, inputCentre, offset, data_region] = ...
    recentre_HPL_dpstack_simple(dat, newCentrePoint, set_bound_value)
% Function to recentre a 'HPL format' dp stack (slices specified in z-direction),
% using a simple (intger value) data shift. Values outside the original data space are
% assigned the value given in set_bound_value.
%
% Usage:
%
% function [circshifted_dat, inputCentre, offset, data_region] = ...
%       recentre_HPL_dpstack_simple(dat, newCentrePoint, set_bound_value)
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

if nargin < 3
    set_bound_value = [];
end

Nin = size(dat,[2 3]);
Nout = Nin;
inputCentre = floor(Nin/2)+1;

newCentrePointUsed = round(newCentrePoint);
offset = inputCentre - newCentrePointUsed;

R1 = max(offset(1)+1,1);
R2 = min(offset(1)+Nin(1),Nout(1));
C1 = max(offset(2)+1,1);
C2 = min(offset(2)+Nin(2),Nout(2));
data_region = [R1, R2, C1, C2];

circshifted_dat = circshift(dat, [0, offset]);

if ~isempty(set_bound_value)
    % set the wrapped over values to set_bound_value:
    if offset(1) < 0
        circshifted_dat(:, (R2 + 1):end,:) = set_bound_value;
    end
    if offset(1) > 0
        circshifted_dat(:, 1:(R1-1),    :) = set_bound_value;
    end
    if offset(2) < 0
        circshifted_dat(:, :, (C2 + 1):end) = set_bound_value;
    end
    if offset(2) > 0
        circshifted_dat(:, :, 1:(C1-1)) = set_bound_value;
    end
end
% Could also do the above by assigning a big zeros or ones block and assigning data region within
% that cube... that might be quicker?


end

