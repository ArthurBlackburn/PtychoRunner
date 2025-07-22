function [included_angle, direction] = IncAngle(alpha, beta, degs)
% Used to find the smallest (included) angle between two angles.
%
% [included_angle, direction] = IncAngle(alpha, beta, degs)
% 
% If the arguments are in degrees, set 'degs' not equal to zero.
% 
% direction: Gives the sense of the included angle from alpha to beta. If
% 0 the direction is anticlockwise.
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

if nargin == 2
    degs = 0;
end

if degs
    alpha = (alpha/180)*pi;
    beta = (beta/180)*pi;
end

alpha_sign = sign(alpha);
beta_sign = sign(beta);

if beta_sign ~= alpha_sign
    if beta_sign<0
        beta = 2*pi + beta;
    else
        beta = -2*pi + beta;
    end
    
end

included_angle = beta - alpha;

if alpha_sign < 0
    if (-pi < included_angle) &&(included_angle < 0)
        direction = 1;
        included_angle = -included_angle;
    else
        direction = 0;
        if included_angle < -pi
            included_angle = included_angle+(2*pi);        
        end
    end
else
    if (0 < included_angle) &&(included_angle < pi)
        direction = 0;
    else
        direction = 1;
        if included_angle > pi
            included_angle = 2*pi - included_angle;
        else
            included_angle = -included_angle;
        end
    end
end

if degs
    included_angle = (included_angle/pi)*180;
end

end


