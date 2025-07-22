function data = loadHPL(filename, path_in_file)
%
% Loads data from HPL (HDF5) format used by Azorus, or other programs, 
% placing datasets and contained JSON format meta data within structs.
%
% Usage:
% data = loadHPL(filename)
% data = loadHPL(filename, path_in_file)
% 
% Parameters
% ----------
%
% filename
%     Name of the file to load data from
% path_in_file : optional
%     Path to the part of the file to load. E.g. for just coords, enter
%     'coords', or for the for the first step of reconstruction 
%     'pi_recipe/r0'
%
% Output data contains from Azorus contains amongst other hopefully self-
% explanatory fields, 'meta' which has fields of note:
%
% data.meta.angular_ps: The angular pixel size (radian) 
% data.meta.affineRotation: The rotation between the diff_patts and coords.
%   This rotation will have already been applied to the coords set, to align
%   dp to recon space. If you want to rotate recon etc back at end then  
%   this rotation in radians will need to be applied.
%
% Note also that:
% 'coords_nm' are the calibrated co-ordinates in nm, where the first column
% is x and the second is y, in a conventional cartesian system ('x across, y is
% up).
% 'coords' are the voltage coords in the scan where the first column is 'y',
% where y increases in value from the top downwards, thus akin to row number 
% in image / matrix, and the second column is the x-scan voltage, akin to 
% column number.
%
% HDF5 import below largely based on code from Pauli Virtanen <pav@iki.fi>,
% as uploaded to Matlab file exchange, described as in the Public Domain 
% with no warranty. Adapted and updated by Arthur Blackburn, University of
% Victoria, 2020, <ablackbu@uvic.ca>

if nargin > 1
  path_parts = regexp(path_in_file, '/', 'split');
else
  path_in_file = '';
  path_parts = [];
end
loc = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
try
  data = load_one(loc, path_parts, path_in_file);
  H5F.close(loc);
catch exc
  H5F.close(loc);
  rethrow(exc);
end

% give the (hopefully) calibrated coords:
if isfield(data,'coords') && isfield(data,'meta')
    % swap x and y, swap sign of y and mult by calibration (nm V^-1)
    data.coords_nm = data.coords(:,2:-1:1) * [[1 0];[0 -1]] * data.meta.scanCalibration*1e9;
end

end

function data=load_one(loc, path_parts, full_path)
% Load a record recursively.
while ~isempty(path_parts) && strcmp(path_parts{1}, '')
  path_parts = path_parts(2:end);
end
data = struct();
num_objs = H5G.get_num_objs(loc);
% 
% Load groups and datasets
%
for j_item=0:num_objs-1
  objtype = H5G.get_objtype_by_idx(loc, j_item);
  objname = H5G.get_objname_by_idx(loc, j_item);
  
  if objtype == 0
    % Group
    name = regexprep(objname, '.*/', '');
  
    if isempty(path_parts) || strcmp(path_parts{1}, name)
      if ~isempty(regexp(name,'^[a-zA-Z].*','once'))
          
        group_loc = H5G.open(loc, name);
        try
          sub_data = load_one(group_loc, path_parts(2:end), full_path);
          H5G.close(group_loc);
        catch exc
          H5G.close(group_loc);
          rethrow(exc);
        end
        if isempty(path_parts)
          data.(name) = sub_data;
        else
          data = sub_data;
          return
        end
      end
    end
   
  elseif objtype == 1
    % Dataset
    name = regexprep(objname, '.*/', '');
  
    if isempty(path_parts) || strcmp(path_parts{1}, name)
      if ~isempty(regexp(name,'^[a-zA-Z].*','once'))
        dataset_loc = H5D.open(loc, name);
        try
            sub_data = H5D.read(dataset_loc, ...
              'H5ML_DEFAULT', 'H5S_ALL','H5S_ALL','H5P_DEFAULT');
            H5D.close(dataset_loc);
        catch exc
            H5D.close(dataset_loc);
            rethrow(exc);
        end
	
        sub_data = fix_data(sub_data, name);
	
        if isempty(path_parts)
          data.(name) = sub_data;
        else
          data = sub_data;
          return
        end
      end
    end
  end
end

% Check that we managed to load something if path walking is in progress
if ~isempty(path_parts)
  error('Path "%s" not found in the HDF5 file', full_path);
end

end

function data=fix_data(data, name)
% Fix some common types of data to more friendly form.
if isstruct(data)
  fields = fieldnames(data);
  if length(fields) == 2 && strcmp(fields{1}, 'r') && strcmp(fields{2}, 'i')
    if isnumeric(data.r) && isnumeric(data.i)
      data = data.r + 1j*data.i;
    end
  end
end

if ~strcmpi(name,'diffPats')
    if isnumeric(data) && ndims(data) > 1
      % permute dimensions
      data = permute(data, fliplr(1:ndims(data)));
    end
end

if ischar(data) && (length(data) > 1)
    try
        if iscolumn(data)
            data = data';
        end
        data = jsondecode(data);
    catch
        disp('INFO: non-JSON format char strings encountered');
    end
end

end