function [exit_wave, nslices] = propogate_object(obj, obj_info, entrance_wave_func)
% Used to propogate an entrance wavefunction through a multi-slice object model.
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
    entrance_wave_func = [];
end

if iscell(obj)
    obj_size = size(obj{end});
    obj = obj{end};
else
    obj_size = size(obj);
end


if length(obj_size) > 2
    nslices = obj_size(end);
else
    nslices = 1;
end

if (isfield(obj_info,'delta_z_m') && ~isfield(obj_info,'delta_z_angstrom')) || ...
    isfield(obj_info,'delta_z_m') && (isfield(obj_info,'delta_z_angstrom') && isempty(obj_info.delta_z_angstrom))
    obj_info.delta_z_angstrom = obj_info.delta_z_m/1e-10;
end

if nslices > 1
    obj_inds = num2cell(ones(1,length(obj_size)));
    obj_inds(1:2) = {':'};
    im0 = squeeze(obj(obj_inds{:}));
    if ~isempty(entrance_wave_func)
        im0_size = size(im0);
        [X, Y] = meshgrid((1:1:im0_size(2))*obj_info.recon_pixel_size, ...
                          -1*((1:1:im0_size(1))*obj_info.recon_pixel_size));
        Z_ent = 0;
        ent_wave = entrance_wave_func(X, Y, Z_ent);
        Z_exit = (nslices - 1) * obj_info.delta_z_angstrom*1e-10;
        unmodified_exit_wave = entrance_wave_func(X, Y, Z_exit);
        im0 = im0.*ent_wave;
    end
    for slc = 2:obj_size(end)
        im1_illum = slice_progate(im0, ...
            obj_info.recon_pixel_size,...
            obj_info.lambda, ...
            obj_info.delta_z_angstrom*1e-10);
        obj_inds{end} = slc;
        im0 = im1_illum .* squeeze(obj(obj_inds{:}));
    end    
    if ~isempty(entrance_wave_func)        
        exit_wave = im0./unmodified_exit_wave;
    else
        exit_wave = im0;
    end
else
    exit_wave = obj;
end


end

