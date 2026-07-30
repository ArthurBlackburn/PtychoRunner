%hdf5_pos loads scan positions from hdf5 files 
%Written by Yi Jiang & Zhen Chen,
% Modified by Arthur Blackburn, UVic, 2021 onwards

function [ p ] = hdf5_pos( p )

for ii = 1:p.numscans
%     positions_real = zeros(0,2);

    if isfield(p.scan,'preloaded_positions_real') && ~isempty(p.scan.preloaded_positions_real)
        % Added by A. Blackburn: PtychoProcessFromTableF_MS.m already parsed
        % the source file (.mat/.hdf5/.h5) and applied trans_coords and
        % rotation (including drift correction and subregion cropping, if
        % any) to produce ps.dps.rotated_posns, so reuse that directly rather
        % than re-reading the file and re-deriving positions here. The
        % trans_coords/scale_in_ps/rotation steps below are intentionally
        % skipped in this branch - they're already baked into the preloaded
        % array. See the matching preload wiring in Runner/ptychoset_mult.m's
        % conv_ptychoset_cSAXS_MSlc and the equivalent short-circuit in
        % +detector/+prep_data/+matlab_aps/load_data.m.
        utils.verbose(2, 'Using pre-loaded scan positions from PtychoProcessFromTableF_MS.m; skipping file read, trans_coords, and rotation (already applied upstream).');
        positions_real = p.scan.preloaded_positions_real;
    else

    switch p.scan.type
        case 'custom'
            if isempty(p.scan.custom_positions_source) %guess the position file name from base path
                pos_file = strcat(p.base_path,sprintf(p.scan.format, p.scan_number(ii)),'/data_position.hdf5');
            else
                pos_file = p.scan.custom_positions_source; % added by Arthur B
            end
            if exist(pos_file,'file')
                [~, ~, crds_EXT] = fileparts(pos_file);
                if strcmpi(crds_EXT,'.hdf5')                
                    ps = h5read(pos_file,'/probe_positions_0');
                    ppX = ps(:,1)*1e-10; % Angstrom to meter --> ppX is read in column 1
                    ppY = ps(:,2)*1e-10; % --> ppY is read in column 2
                elseif strcmpi(crds_EXT,'.mat')
                    % ... a mat file also be read as hdf5 if is is saved as v7.3...
                    ps = load(pos_file);
                    ppX = ps.probe_positions_0(:,1)*1e-10; % Angstrom to metre --> ppX is read in column 1
                    ppY = ps.probe_positions_0(:,2)*1e-10; % --> ppY is read in column 2
                elseif strcmpi(crds_EXT,'.h5')
                    % py4DSTEM/EMD file produced by ptyzer's azohp_to_py4d.py
                    % (D:\Data\Repos\ptyzer). Added by A. Blackburn, to allow
                    % positions to be read straight from the same combined .h5
                    % file the diffraction data is loaded from (see the matching
                    % case '.h5' in +detector/+prep_data/+matlab_aps/load_data.m).
                    % Positions are read from the already-flattened
                    % scan_coords_pointlist group (raveled in numpy order='C',
                    % i.e. the same flat scan-position order that
                    % load_data.m's '.h5' case reproduces for the diffraction
                    % stack) rather than the 2D scan_coords grid, so no
                    % row/column ravel-order ambiguity has to be resolved here.
                    % Units are nm (already-calibrated real-space coordinates),
                    % not Angstrom/volts, and x/y are not pre-swapped the way
                    % ppX/ppY are for the other two cases below.
                    ppX = h5read(pos_file, '/py4DSTEM_root/datacube/scan_coords_pointlist/x');
                    ppY = h5read(pos_file, '/py4DSTEM_root/datacube/scan_coords_pointlist/y');
                    ppX = ppX(:)*1e-9; % nm to metre
                    ppY = ppY(:)*1e-9; % nm to metre
                    % azohp_to_py4d.py's scan_coords are already real-space,
                    % Cartesian-oriented (x,y) nm coordinates (unlike the raw
                    % voltage 'coords' field the other two cases ultimately
                    % derive from), so trans_coords = 3 (identity: [ppX, ppY])
                    % is the most likely starting point - but verify visually
                    % against a known scan direction before trusting it, same
                    % as for any other position source here.
                else
                    error('unrecognized postiion file type')
                end
                positions_real = zeros(length(ppX),2); 

                positions_real(:,1) = -ppY(:); % Note column 1 is -1*read in column 2
                positions_real(:,2) = -ppX(:); % Note column 2 is -1*read in column 1 , 
                % This would be the same as reading do h5read then doing transcoords case 6 
                % in PtychoProcessFromTableF (Comment by AB).                
            else
            	error('Could not find function or data file %s', pos_file);
            end
        otherwise
            error('Unknown scan type %s.', p.scan.type);
    end

    if isfield(p.scan,'trans_coords') && ~isempty(p.scan.trans_coords)
        switch(p.scan.trans_coords)
            case 1
                positions_real = [-ppY,  ppX];
            case 2
                positions_real = [ ppY, ppX]; %
            case 3
                positions_real = [ ppX,  ppY]; %
            case 4 
                positions_real = [ ppX, -ppY]; % 
            case 5 
                positions_real = [ ppY, -ppX]; % if ppX and ppY are cartesian coords this seems to work. 
                % Thus it would appear that positions_real coord system is rotated by 90 from my input cartesian
            case 6
                positions_real = [-ppY, -ppX]; % the default for files from Azorus... at the moment..
            case 7 
                positions_real = [-ppX, ppY]; % 
            case 8 
                positions_real = [-ppX, -ppY]; % 
            otherwise
                error('Invalid trans_coords');
        end
        utils.verbose(2, 'Transposed hdf5 read coordinates in hdf5_pos case: %05d.', p.scan.trans_coords)
    end
    
    if isfield(p.scan,'scale_in_ps') && ~isempty(p.scan.scale_in_ps) && ~isnan(p.scan.scale_in_ps)
        % if scale in ps is used, assume that we have been passed voltage
        % coordinates. Above assume Angstroms (10^-10m) were coming in
        % scale_in_ps and nm_per_V multiliply together with above read in coords to give co-ords in metres...
        if length(p.scan.scale_in_ps) == 1
            scale_xy = p.scan.scale_in_ps*[1 1];
        end
        positions_real = scale_xy .* positions_real ;
    end
    
    % After all this, positions real ends up being [Y, -X]
    
    if isfield(p.scan,'rotation') && ~isempty(p.scan.rotation)
        
           rot_angle = (p.scan.rotation /180)*pi;

           R = [cos(rot_angle) -sin(rot_angle);...
                sin(rot_angle)  cos(rot_angle)];
            
           rotated_posns = R*(positions_real.');
           positions_real = rotated_posns.';
    end

    end

    p.numpts(ii) = size(positions_real,1);
    p.positions_real = [p.positions_real ; positions_real]; %append position - 
    % comment AB: 
    % this reveals the nature of num scans: i.e. could have multiple files for
    % postions for scans, which are all loaded up. This is not usually the
    % case in electron ptychography. Though might eventually need something like
    % this to deal with iterations and dose fractionation etc.

    
end
    
end

