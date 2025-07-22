function z = circ(x, y, D)
% A very simple function for creating a circular 'mask' based on x, y and Diameter of circle
% 
% function z = circ(x, y, D)
% 
% z = double(((x.^2+y.^2)<D/2) except if sqrt(x.^2+y.^2) == D/2, which gives z = 0.5.
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

r2 = x.^2 + y.^2;
z = double(r2 < ((D^2)/4));
z(r2 == ((D^2)/4)) = 0.5;