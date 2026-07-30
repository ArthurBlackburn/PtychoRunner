%LOAD_DATA prepare filenames and load data
% ** p      p structure
%
% returns:
% ++ p      p structure
%
% see also: detector.prep_data.matlab_ps.prepare_data
% written by YJ


function [ p ] = load_data( p )
import utils.find_files
import io.image_read

det.params = p.detectors(p.scanID).params;
detStorage = p.detectors(p.scanID).detStorage;
%% prepare filenames
if isempty(detStorage.files)
    [p] = det.params.get_filename(p);
    if numel(detStorage.files) == 0
        error('Did not find any files.')
    end
end

utils.verbose(2, 'Loading raw data of scan %05d.', p.scan_number(p.scanID))

files = detStorage.files;
[~,~,ext]=fileparts(files{p.scanID});
ismask = false;
if isfield(p,'preload') && isfield(p.preload,'diffraction') && ~isempty(p.preload.diffraction)
    % Added by A. Blackburn: PtychoProcessFromTableF_MS.m already parsed the
    % source file (.mat/.hdf5/.h5) and applied subregion cropping / dose
    % adjustment, so reuse that array directly rather than re-reading and
    % re-deriving it from disk here. See the matching preload wiring in
    % Runner/ptychoset_mult.m's conv_ptychoset_cSAXS_MSlc and the equivalent
    % short-circuit in +scans/+positions/hdf5_pos.m.
    utils.verbose(2, 'Using pre-loaded diffraction data from PtychoProcessFromTableF_MS.m; skipping file read.');
    data = p.preload.diffraction;
else
switch lower(ext)
    case '.hdf5'
        data = h5read(files{p.scanID},'/dp');
    case '.mat'  % add .mat format by Zhen Chen
        % added by ABlackburn to allow mask file to be embedeed in data if required.
        data=load(files{p.scanID});
        if isfield(data,'mask')
            mask = ~data.mask;
            ismask = true;
        end
        % data=load(files{p.scanID},'dp');
        data=double(data.dp);
    case '.h5'
        % py4DSTEM/EMD file produced by ptyzer's azohp_to_py4d.py (D:\Data\Repos\ptyzer).
        % Added by A. Blackburn, so that Azorus data converted straight to py4DSTEM's
        % format can be reconstructed without first exporting to the separate .mat +
        % position .hdf5 pair PrepareData.m produces. See +scans/+positions/hdf5_pos.m
        % for the matching scan-position loader (case '.h5' there).
        h5file = files{p.scanID};
        h5root = '/py4DSTEM_root/datacube';
        % grid_shape = [Nrows, Ncols], the Azorus scan mesh shape. Scan positions were
        % raveled into scan_coords_pointlist (see hdf5_pos.m) in the same order
        % (numpy order='C', i.e. Ncols/columns fastest) as the original Azorus
        % acquisition, so that order must be reproduced here for the diffraction
        % stack to stay aligned, point for point, with the loaded positions.
        grid_shape = h5read(h5file, [h5root '/metadatabundle/general_meta/meshParams/shape']);
        data4d = h5read(h5file, [h5root '/data']);
        % data4d was written by h5py/emdfile with numpy shape
        % (Rx, Ry, Qx, Qy) = grid_shape + [detector pixel dims]. MATLAB's HDF5
        % readers (h5read, and the low-level H5D.read used in loadHPL.m's
        % fix_data()) return array dimensions in REVERSED order relative to how
        % they were written, so data4d actually comes back sized [Qy, Qx, Ry, Rx].
        % Collapsing the trailing two ([Ry, Rx]) dimensions with a plain reshape
        % (no permute needed) reproduces the same flat scan-position order as
        % scan_coords_pointlist / the original Azorus acquisition, since MATLAB's
        % native column-major flattening of a fully dimension-reversed array walks
        % the same underlying bytes, in the same order, as numpy's default
        % (order='C') ravel of the original array.
        sz = size(data4d);
        n_scan = prod(double(grid_shape));
        if numel(sz) ~= 4 || prod(double(sz(3:4))) ~= n_scan
            error(['Unexpected shape for %s: got %s, expected the last two ', ...
                   'dimensions to multiply out to %d scan positions (grid_shape = [%d %d]). ', ...
                   'The assumed reversed-axis-order convention (see comments above) may not ', ...
                   'hold for this file.'], ...
                   [h5root '/data'], mat2str(sz), n_scan, grid_shape(1), grid_shape(2));
        end
        data = reshape(data4d, sz(1), sz(2), n_scan);
        clear data4d
        % data is now [Qy, Qx, n_scan]: the two detector-pixel axes come back
        % transposed relative to Azorus's native axis order, since py4DSTEM's
        % DataCube labels its detector axes 'Qx'/'Qy' purely positionally, not
        % semantically. Verify this (e.g. against a virtual bright-field image, or
        % against the mean diffraction pattern already stored at
        % [h5root '/dp_mean/data'] in the same file) and, if it looks transposed,
        % set det.params.orientation(1) = true for this detector - that transpose
        % is then applied generically by the orientation-handling code below.
    otherwise
        error('data format not supported!');
end
end
data = squeeze(data);
utils.verbose(2, 'Loaded data from: %s', files{p.scanID});

if det.params.orientation(1)
    utils.verbose(2, 'Transposing diffraction patterns')
    data = permute(data, [2 1 3]);
    if ismask
        mask = permute(mask, [2 1 3]);
    end
end

if det.params.orientation(2) && det.params.orientation(3)
    utils.verbose(2, 'Flipping diffraction patterns')

    data = rot90(data,2);  % merge the fliplr and flipup operations
    if ismask
        mask = rot90(mask,2);
    end
else
    if det.params.orientation(2)
        data = fliplr(data);
        if ismask
            mask = fliplr(mask);
        end
    end
    if det.params.orientation(3)
        data = flipud(data);
        if ismask
            mask = flipud(mask);
        end
    end
end

%{
%% select data loading routine
current_version = version;
ver_str = strsplit(current_version, '.');
ver_num = str2double([ver_str{1} '.' ver_str{2} ver_str{3}]);
if ver_num >= 9.4 
    c_reader = true;
else
    utils.verbose(3, 'Fast data reader is not available for Matlab version %d. Switching to image_read.', ver_num)
    c_reader = false;
end
if c_reader && ~ismember(det.params.file_extension, {'h5', 'cbf', 'tiff', 'tif'})
    utils.verbose(3, 'Fast data reader is not available for selected file format %s. Switching to image_read.', det.params.file_extension)
    c_reader = false;
end


%% load data
if c_reader
    % use fast data reader
    if strcmpi(det.params.file_extension, 'tif')
        det.params.file_extension = 'tiff'; 
    end
    if strcmpi(det.params.file_extension, 'tiff')
        % keep results consistent with matlab's imread and io.image_read for tiff files 
        det.params.orientation(1) = ~det.params.orientation(1); 
    end
    % convert PtychoShelves center to raw data center
    arg.ctr = detStorage.ctr-1;
    sz = det.params.geometry.sz;

    if det.params.orientation(3)
        arg.ctr(1) = round(sz(1) - arg.ctr(1));
    end
    if det.params.orientation(2)
        arg.ctr(2) = round(sz(2) - arg.ctr(2));
    end  
    if det.params.orientation(1)
        arg.ctr = fliplr(arg.ctr);
    end
    if iscolumn(arg.ctr)
        arg.ctr = arg.ctr';
    end
    
    % flip XY
    arg.ctr = fliplr(arg.ctr);
    
    % create structure for c_reader
    arg.data_path = files;
    arg.nthreads = p.io.data_nthreads;
    arg.precision = p.io.data_precision;
    arg.extension = det.params.file_extension;
    arg.asize = detStorage.read_size;
    arg.data_location = detStorage.h5_group;
    assert(all(arg.ctr>0),'Raw data center position has to be positive.')
    % load data and permute
    utils.verbose(2, 'Loading raw data of scan S%05d.', p.scan_number(p.scanID))
    data = io.read_measurement(arg);
    data = squeeze(data);
    size(data)

    if det.params.orientation(1)
        data = permute(data, [2 1 3]);
    end
    if det.params.orientation(2) && det.params.orientation(3)
        data = rot90(data,2);  % merge the fliplr and flipup operations 
    else
        if det.params.orientation(2)
            data = fliplr(data);
        end
        if det.params.orientation(3)
            data = flipud(data);
        end
    end
    
    
    
else
    if isempty(detStorage.h5_group)
        if numel(files)==1
            % one directory contains single file
            dataaux = image_read(files{1}, det.params.image_read_extraargs);
            data = dataaux.data;
        else
            % multiple files per directory
            dataaux = image_read(files, det.params.image_read_extraargs);
            data = dataaux.data;
        end
    else
        data = zeros([p.asize numel(detStorage.h5_group)]);
        if numel(files)==1
            if numel(detStorage.h5_group)==1
                numel(detStorage.h5_group)
                % one hdf5 file; one group
                dataaux = image_read(files{1}, det.params.image_read_extraargs{:}, 'H5Location', detStorage.h5_group{1});
                data = dataaux.data;
            else
                % one hdf5 file; multiple groups
                for ii=1:numel(detStorage.h5_group)
                    dataaux = image_read(files{1}, det.params.image_read_extraargs{:}, 'H5Location', detStorage.h5_group{ii});
                    data(:,:,ii) = dataaux.data;
                end
            end
        else
            if numel(detStorage.h5_group)==1
                % multiple files, single H5 group
                dataaux = image_read(files, det.params.image_read_extraargs{:}, 'H5Location', detStorage.h5_group{1});
                data = dataaux.data;
            else
                % multiple files, but different H5 group
                for ii=1:numel(detStorage.h5_group)
                    dataaux = image_read(files{ii}, det.params.image_read_extraargs{:}, 'H5Location', detStorage.h5_group{ii});
                    data(:,:,ii) = dataaux.data;
                end
            end
        end
    end
end
%}

detStorage.data = data;

if ismask
    detStorage.mask = mask;
end

end
