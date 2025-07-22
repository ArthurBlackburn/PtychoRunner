function keyboard_m(verbose_level, keyboard_thresh)
% Executes a keyboard pause if the verbose level is greater than 
% keyboard_thresh
% Arthur Blackburn, UVic, 2021
%
% Can subsitute keyboard commnand, to allow uninterupted flow of programs 
% even if error conditions occur while in a debug mode. 
% 
% e.g. in PtychoShelves
% search / replace:
% keyboard / utils.keyboard_m(utils.verbose())

if nargin < 2
    keyboard_thresh = 4;
end

if verbose_level > keyboard_thresh
    keyboard;
end

end

