% INITIALIZE_SOLVER initialize GPU ptycho reconstruction, generate cache values, fftshift data, etc 
% 
% [self, cache] = initialize_solver(self,par) 
% 
% ** self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ** par       structure containing parameters for the engines 
%
% returns: 
% ++ self      structure containing inputs: e.g. current reconstruction results, data, mask, positions, pixel size, ..
% ++ cache     structure with precalculated values to avoid unnecessary overhead

% Academic License Agreement
% 
% Source Code
% 
% Introduction 
% �	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
% 
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE�s responsibility to ensure its proper use and the correctness of the results.�
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
% 
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379�382 (2008). 
%   (doi: 10.1126/science.1158573),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68�71 (2013). (doi: 10.1038/nature11806),
% for LSQ-ML method 
% M. Odstrcil, A. Menzel, M.G. Sicairos,  Iterative least-squares solver for generalized maximum-likelihood ptychography, Optics Express, 2018
% for OPRP method 
%  M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  "Ptychographic coherent diffractive imaging with orthogonal probe relaxation." Optics express 24.8 (2016): 8360-8369
% and/or for multislice:
% E. H. R. Tsai, I. Usov, A. Diaz, A. Menzel, and M. Guizar-Sicairos, X-ray ptychography with extended depth of field, Opt. Express 24, 29089�29108 (2016). 
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       � All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Z�rich, Switzerland. 
% 
%   

function [self, cache] = init_solver(self,par)
        
    import engines.GPU.shared.*
    import math.*
    import utils.*
    import plotting.*
    import engines.GPU.GPU_wrapper.*

    par.Nscans = length(self.reconstruct_ind);
    cache.skip_ind = setdiff(1:self.Npos,[self.reconstruct_ind{:}]); %  wrong datasets  to skip 
        
    if ~any(self.probe_support(:))
        self.probe_support = []; 
    end
    %% avoid probe to be larger than a certain oversampling !!!! 
    if isempty(self.probe_support) 
        par.probe_backpropagate = 0;
    end

    if  ~isempty(self.background) && any(self.background(:) > 0)
        Background = self.background;
    elseif  par.background_detection
        Background = 0; 
    else
        Background = [];  %  array of background light 
    end
    
    Noise = [];
    %% prepare data / noise / mask 
    if par.relax_noise &&  ~isempty(self.noise) &&  strcmp(par.likelihood, 'L1')
        Noise = self.noise;
%       Noise = (sqrt(posit(self.diffraction + Noise)) - sqrt(posit(self.diffraction - Noise)))/2;
% In the above ^^^
% posit is not a built-in matlab function. I think the intention here might be to only let
% posit's argument pass only if it is _posit_ive??, i.e. if -ve set to zero, so we have 
% no complex numbers arising from the below
        Noise = (sqrt(max(self.diffraction + Noise,0)) - sqrt(max(self.diffraction - Noise,0)))/2;
        Noise(self.diffraction == 0) = 1;
        disp('Using measured noise')
        Noise = max(0.5, Noise); 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% PREPARE MASK AND DATA %%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %% prepare mask , note that bool in matlab has size of uint8 !!
    cache.mask_indices = [];
    if any(self.mask(:))
        Mask = [];
        % single mask 
        if all(all(mean(self.mask,3) == self.mask(:,:,1)))
            Mask = self.mask(:,:,1);
        else
            % mask for each scan
            for ll = 1:par.Nscans
                ind = self.reconstruct_ind{ll};
                %if there is only one repeated mask over whole scan 
                if all(all(all(bsxfun(@eq, self.mask(:,:,ind), self.mask(:,:,ind(1))))))
                    Mask(:,:,ll) = self.mask(:,:,ind(1)); 
                end
                cache.mask_indices(ind) = ll;
            end
        end
        if isempty(Mask)
            % mask for each position 
            Mask = self.mask;  % otherwise just store original
            cache.mask_indices(ind) = 1:self.Npos;
        end
        % important to save memory 
        if all(Mask(:) == 1 | Mask(:) == 0)
            Mask = logical(Mask );
        else
            Mask = uint8(Mask*255);  % if there are nonlogical values in mask, store them as uint8 to save memory 
        end
    else
        Mask = [];
    end
    
    % Arthur Blackburn, University of Victoria,
    % Added code to deal with creating masks
    if ~isempty(par.p.border_mask) && par.p.border_mask(1) > 0
        MaskBM = zeros(par.p.asize,'logical');
        if length(par.p.border_mask) == 1
            % Apply a simple boundary mask of width xx to the data.
            mw = par.p.border_mask;
            MaskBM([1:mw,(end-(mw-1)):end], :) = 1;
            MaskBM(:, [1:mw,(end-(mw-1)):end]) = 1;
            MaskBM = fftshift_2D(MaskBM);
        elseif length(par.p.border_mask) == 3
            % assume it's a cicular mask we are using here, for example to
            % mask out the region of the HAADF detector in the SU9000 data.
            if par.p.border_mask(1) == 0
                r_cen = round(par.p.asize(1)/2)+1;
            else
                r_cen = par.p.border_mask(1);
            end
            if par.p.border_mask(2) == 0
                c_cen = round(par.p.asize(2)/2)+1;
            else
                c_cen = par.p.border_mask(2);
            end
            radius = par.p.border_mask(3);
            x_pts = (-r_cen-1):1:(par.p.asize(1)-r_cen-2);
            y_pts = (-c_cen-1):1:(par.p.asize(2)-c_cen-2);
            [xx, yy] = meshgrid(x_pts, y_pts);
            MaskBM = (xx.^2 + yy.^2) > radius^2;
            MaskBM = fftshift_2D(MaskBM);
            clear xx yy
        else
            if ischar(par.p.border_mask) && strcmpi(par.p.border_mask, 'imdef')
                % Define the mask by looking at values which are
                % exactly zero in the sum of the all the diffraction images.
                dp_sum = squeeze(sum(self.diffraction,3));
                Mask = dp_sum == 0;                
            else
                % Assume that a 4 element array where
                % corner points are: upper row, lower row, left col, right col.
                MaskA = MaskBM; MaskB = MaskBM;
                MaskA(   par.p.border_mask(1):par.p.border_mask(2), :) = 1;
                MaskB(:, par.p.border_mask(3):par.p.border_mask(4))    = 1;
                MaskBM = ~(MaskA & MaskB);
                if length(par.p.border_mask) == 5
                    % Then the 5th element is used to specify pincusion distortion...
                    % However, 'imdef' above would also work well in most experimental cases 
                    % for pincusion distortion.
                    lens_options = {'borderType','fix',...
                    'interpMethod', 'linear',...
                    'padMethod','fill',...
                    'padValue',1,...
                    'fType',4,...
                    'intensCorr',false};
                    MaskBM = lensdistort(single(MaskBM), par.p.border_mask(5), lens_options{:});
                    MaskBM = logical(MaskBM);
                end
                MaskBM = fftshift_2D(MaskBM);
                clear MaskA MaskB
            end
        end
        % self.mask and Mask might have been scaled in other functions:
        if any(par.p.asize_presolve ~= par.p.asize)
            fill_value = 1;  % fill value in case of ptychographic "super resolution"
            MaskBM = fftshift_2D(MaskBM);
            MaskBM = crop_pad(MaskBM, self.Np_p, fill_value);
            MaskBM = ifftshift_2D(MaskBM);                        
        end                
        if ~isempty(Mask)
            Mask = Mask | MaskBM;
        else
            Mask = MaskBM;
        end
        clear MaskBM;
%         figure; imagesc(Mask); axis image;
    end
    
    if ~isempty(par.p.cross_mask) && par.p.cross_mask(1) > 0
        MaskBM = zeros(par.p.asize,'logical');
        MaskA = MaskBM; MaskB = MaskBM;
        MaskA(   par.p.cross_mask(1):(par.p.cross_mask(1) + par.p.cross_mask(3)), :) = 1;
        MaskB(:, par.p.cross_mask(2):(par.p.cross_mask(2) + par.p.cross_mask(3)))    = 1;
        MaskBM = MaskA & MaskB;
        MaskBM = fftshift_2D(MaskBM);
        clear MaskA MaskB
        if ~isempty(Mask)
            Mask = Mask | MaskBM;
        else
            Mask = MaskBM;
        end
        clear MaskBM;
    end

    
    
    %% prepare diffraction data 
    Diffraction = self.diffraction;  % diffraction is intensity, not amplitude, comment by Zhen Chen


    if par.upsampling_data_factor
        % downsample the data down to original size to save memory 
        Diffraction = utils.binning_2D(Diffraction, 2^par.upsampling_data_factor) * (2^(2*par.upsampling_data_factor)); 
        if ~isempty(Mask)
            Mask = utils.binning_2D(Mask, 2^par.upsampling_data_factor) == 1; 
        end
    end
    
    Diffraction = single(max(0,Diffraction));

    
    if ~isempty(Mask)
        if size(Mask,3) == par.Nscans && par.Nscans  > 1
            for ll = 1:par.Nscans
                ind = self.reconstruct_ind{ll};
                Diffraction(:,:,ind) = Diffraction(:,:,ind) .* ~Mask(:,:,ll); 
            end
        else
            Diffraction = Diffraction .* ~Mask;
        end
    end
    
      
    
    if  ~isinf(self.z_distance(end)) %  && mod(Ninf,2)~=0
        % assume inputs already fftshifted, but in case of nearfield
        % fftshift it back for the ASM propagator 
        Noise = fftshift_2D(Noise);
        Diffraction = fftshift_2D(Diffraction);
        Mask = fftshift_2D(Mask);
    end

    
    %%%% compress data if requested %%%%%%
    if par.compress_data
        DATA_MAX = quantile(max2(abs(Diffraction)), 1-1e-2); 
        C_factor_0 = 2;  % compression factor  >=2 seems to be safe, >=4 is pratically lossless 
        if par.compress_data == 1 || DATA_MAX < 2^(2*8) / C_factor_0^2
            Diffraction = sqrt(single(Diffraction));
            if  DATA_MAX < 2^(2*8) / C_factor_0^2
                % simple sqrt compression to 8 bits 
                verbose(1, 'Online data compression to 8-bits')
                Diffraction = uint8(C_factor_0*Diffraction);
                cache.C_factor = C_factor_0;
            elseif DATA_MAX < 2^(2*16) / 16^2
                % failsafe option: sqrt compression to 16 bits 
                verbose(1, 'Online data compression to 16-bits')
                cache.C_factor = 16;  % use compression factor 16, to be super safe just because we have space
                Diffraction = uint16(cache.C_factor*Diffraction);
            else
                error('Online compression will fail')
            end
        elseif par.compress_data == 2
            % SVD subtraction compression to 8 bits (failsafe is compression to 16bits)
            % additionally remove some SVD modes 
            Diffraction = sqrt(single(Diffraction));
            Nmodes = par.Nscans;
            [U,S,V] = fsvd(reshape(Diffraction,prod(self.Np_p),[]), Nmodes);
            ind_relevant = diag(S).^2/sum(diag(S).^2) > 1e-2; % more than 1% of power
            cache.US_diffraction = (U(:,ind_relevant)*S(ind_relevant,ind_relevant));
            cache.V_diffraction = V(:,ind_relevant);
            svd_Diffraction = round(reshape(cache.US_diffraction*cache.V_diffraction',[self.Np_p, self.Npos]));

            %% compress 
            cDiffraction = single(Diffraction) - svd_Diffraction;
            % reestimate optimal compression factor to keep values < 128 
            C_factor = min(C_factor_0, 128/quantile(max2(abs(cDiffraction)), 1-1e-2));
            if C_factor > 3
                verbose(1, 'Online data compression to 8-bits + SVD')
                cache.C_factor = C_factor;
                cache.US_diffraction =  cache.US_diffraction;
                cache.V_diffraction = cache.V_diffraction*C_factor;
                Diffraction = int8(cDiffraction*C_factor);
            elseif DATA_MAX < 2^(2*16) / 16^2
                % sqrt compression to 16 bits 
                %warning(sprintf('Too high online compression of data, it may cause problems\n Compression factor is %2.2f but should be >= 2\n Switching from 8 to 16bits',C_factor))
                verbose(1, 'Online data compression to 16-bits')
                C_factor = 16;
                cache.C_factor = C_factor;
                Diffraction = uint16(C_factor*Diffraction); 
            else
                error('Online compression will fail')
            end
            
            clear svd_Diffraction cDiffraction 

        else
            error('Unimplented level of compression')
        end
    else
       % precalculate sqrt from the data, store as singles 
        Diffraction = sqrt(single(max(0,Diffraction))); % diffraction is amplitude, comment by Zhen Chen
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% write back the data arrays 
    self.diffraction = Diffraction; % diffraction is amplitude, comment by Zhen Chen
    self.mask = Mask; 
    self.noise = Noise;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% PREPARE GEOMETRY, PROPAGATION, MODES%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % precalculate ASM factor for propagation distance recovery 
    [ASM_difference] = near_field_evolution_gradient(self.Np_p, self.lambda,  self.pixel_size .*self.Np_p );
    cache.ASM_difference = fftshift(ASM_difference);
    
    
    % custom propagator to account for tilted plane sample
    if any(par.p.sample_rotation_angles(1:2)) && check_option(par.p, 'apply_tilted_plane_correction', 'propagation') 
        % get propagators to the tilted plane 
        [tilted_plane_propagate_fwd, tilted_plane_propagate_back] = ...
            get_tilted_plane_propagators(Garray(self.probe{1}), ...
                                        [par.p.sample_rotation_angles(1:2),0],...
                                        self.lambda, self.pixel_size); 
    else
        tilted_plane_propagate_fwd = [];  tilted_plane_propagate_back = [];
    end
    
  
    if ~iscell(self.affine_matrix)
        self.affine_matrix = {self.affine_matrix}; 
    end
    
%     modes = cell(max(par.Nmodes, par.Nlayers),1);
    modes = cell(par.Nlayers,1); % corrected by Zhen Chen
    
%     for i = 1:max(par.Nmodes, par.Nlayers) 
    for i = 1:par.Nlayers % changed by Zhen Chen, Nmodes is the # Nscans.

        verbose(2,'Creating new modes files ')

        modes{i}.lambda =  self.lambda;
        
        
        
        % decompose affine matrix into scale, asymmetry, rotation, shear 
        %affine = scale*[1+asym/2,0; 0,1-asym/2]*[cosd(rot), sind(rot); -sind(rot), cosd(rot)] * [1,0;tand(shear),1];
        
        affine_matrix = self.affine_matrix{min(i,end)};
        [scale, asymmetry, rotation, shear]  = decompose_affine_matrix(affine_matrix);
        
        % store initial geometry parameters 
        modes{i}.scales =    repmat(scale, 1,par.Nscans); 
        modes{i}.asymmetry =  repmat(asymmetry, 1,par.Nscans); 
        modes{i}.shear =   repmat(shear, 1,par.Nscans); 
        modes{i}.rotation =  repmat(rotation, 1,par.Nscans); 
        modes{i}.affine_matrix =  repmat(affine_matrix, 1,1,par.Nscans); 
        modes{i}.shift_scans =  zeros(2, par.Nscans); 
        modes{i}.probe_scale_upd =  0; 
        modes{i}.probe_rotation = ones(1,par.Nscans) * par.sample_rotation_angles(3); % one rotation per scan 
        if par.mirror_objects
            modes{i}.probe_rotation = modes{i}.probe_rotation .* [1,-1];  % flip the coordinates for mirrored object (ie 0 vs 180deg rotation)
        end
        modes{i}.probe_rotation_all = zeros(self.Npos,1);
        for jj = 1:par.Nscans
            modes{i}.probe_rotation_all(self.reconstruct_ind{jj}) = modes{i}.probe_rotation(jj); % one rotation per scan 
        end
        
        distance = self.z_distance(min(end,i)); 
        if ~isinf(distance)
            verbose(2, 'Layer %i distance %g um   ', i, distance*1e6 )
        end
        modes{i}.distances = distance;
            
        if is_used(par, 'fly_scan') && (~isfield(modes{i}, 'probe_positions') || isempty(modes{i}.probe_positions) )
            %% get positions for fly scans 
             self = prepare_flyscan_positions(self, par); 
        else
           %% get positions for normal tomo 
            try  %  try to reuse the positions of there are saved 
                modes{i}.probe_positions = self.modes{i}.probe_positions;
                verbose(2,'Using saved positions')
            catch
                if (modes{i}.scales(end) == modes{1}.scales(end)) && ~isempty(self.probe_positions)
                    modes{i}.probe_positions = self.probe_positions;
                else
                    verbose(2,'Using original positions')
                    modes{i}.probe_positions = (affine_matrix*self.probe_positions_0')';
                end
            end
            try
                modes{i}.probe_positions_0 = self.modes{i}.probe_positions_0;
            catch
                modes{i}.probe_positions_0 = self.probe_positions_0;
            end

        end
        modes{i}.probe_positions_update = { zeros(size(modes{i}.probe_positions)) };
        modes{i}.probe_positions_all = {modes{i}.probe_positions};
        modes{i}.probe_positions_weight = zeros(self.Npos, 1);   
        if isfield(self, 'probe_fourier_shift')&& ~isempty(self.probe_fourier_shift)  && i == 1
            modes{i}.probe_fourier_shift = self.probe_fourier_shift; 
        else
            modes{i}.probe_fourier_shift =  zeros(self.Npos,2);
        end


        if ~isempty(self.probe_support) && i <= par.Nrec
            modes{i}.probe_support = self.probe_support; 
            if i == 1
                verbose(2,'Using real-space probe support')
            end
        else 
            modes{i}.probe_support = [];
        end
        
        if ~isempty(self.probe_support_fft) && i <= par.Nrec && ~check_option(par.p,'probe_support_apt_radius')
            modes{i}.probe_support_fft = fftshift(self.probe_support_fft); 
            if i == 1
                verbose(2,'Using far-field probe support')
            end                      
        elseif check_option(par.p,'probe_support_apt_radius') % not shift for TEM aperture, by Zhen Chen
        % Zhen Chen did not do a fftshift, [but I - A.B. at UVic - am the way I've defined things 
        % (See load_from_p L255)]
            modes{i}.probe_support_fft = fftshift(self.probe_support_fft);
            if i == 1
                verbose(2,'Applying constraint on fft of probe defined by probe_support_apt_radius')
            end 
        else 
            modes{i}.probe_support_fft = [];
        end
        
                                
        F = mean( self.pixel_size)^2 .* mean(self.Np_p) /   (modes{i}.lambda * modes{i}.distances); 
        if F ~= 0
            verbose(3,'Nearfield propagation: Fresnel number/Npix %3.3g', F)
        end
        scale = modes{i}.scales(end);
        modes{i}.ASM_factor = [] ;
        modes{i}.cASM_factor = [] ;

        if ~isinf(modes{i}.distances(end)) % Forward Fresnel propagator in k-space, ASM, commented by Zhen Chen 
            %% near field factor            
            ASM =  exp( modes{i}.distances(end)* cache.ASM_difference);
            modes{i}.ASM_factor = ASM;
            modes{i}.cASM_factor = conj(ASM);
        end
        
        %% far field factor 
        modes{i}.FAR_factor = [];
        modes{i}.cFAR_factor = conj(modes{i}.FAR_factor);
  
        if isinf( par.probe_backpropagate)
             modes{i}.support_fwd_propagation_factor = inf; 
             modes{i}.support_back_propagation_factor = -inf;
        elseif par.probe_backpropagate ~= 0
            [~, modes{i}.support_fwd_propagation_factor] = utils.prop_free_nf( self.probe{1}(:,:,1), par.probe_backpropagate,...
                modes{i}.lambda,  self.pixel_size ./ scale );
  
            modes{i}.support_fwd_propagation_factor = fftshift( modes{i}.support_fwd_propagation_factor );
            modes{i}.support_back_propagation_factor = conj(modes{i}.support_fwd_propagation_factor);
        else
            modes{i}.support_fwd_propagation_factor = [];
            modes{i}.support_back_propagation_factor = [];
        end
        
        modes{i}.tilted_plane_propagate_fwd = tilted_plane_propagate_fwd;
        modes{i}.tilted_plane_propagate_back = tilted_plane_propagate_back;
 
    end
    
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% PREPARE PROBES, INCOHERENT MODES %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    probe_0 = mean(self.probe{1},3);
    probe = cell(par.probe_modes,1);  % force a column cell by Zhen Chen
    for i = 1:par.probe_modes
        try
            probe{i} = self.probe{i};
                %% test if the probe size ok for the variable probe settings etc 
            assert( size(probe{i},4) == 1+par.variable_probe_modes || ...
                ~(par.variable_probe) || i > 1)
            assert((size(probe{i},3) ==1 || par.variable_probe) || ...
                   (size(probe{i},3) == par.Nscans && ~par.share_probe)  )  % no variable prob extension and multiple probes used 
            assert(size(probe{i},3) == par.Nscans || par.share_probe || par.variable_probe, 'Wrong probe size for not shared probe option')
        catch
            if i <= par.Nrec || is_used(par, 'fly_scan')
                verbose(2, 'Creating probe')

                if ~par.share_probe && size(probe{i},3) == 1
                     % dont share probe between scans 
                     probe{i} = repmat(probe_0,[1,1,par.Nscans]);
                end
                if  (par.variable_probe && par.variable_probe_modes > 0) && i == 1
                    verbose(2,'Creating variable probe ')
                    probe{i}(:,:,:,2:1+par.variable_probe_modes) = ...
                        randn([self.Np_p, size(probe{i},3), par.variable_probe_modes])+randn([self.Np_p,size(probe{i},3), par.variable_probe_modes])*1i;
                    continue 
                end
            end
            if length(probe) < i   % none of above 
               % simply create slightly shifted modes in fourier domain, it is useful for
               % inital guess of incoherent modes after orthogonalization 
                step = median(diff(self.probe_positions_0)); 
                probe{i} = 0.01*fftshift(imshift_fft(fftshift(probe_0), randn, randn, false));
            end
            
            % fill the unreconstructed positions if the OPRP method is used
            if par.variable_probe && is_method(par, 'PIE') && i ==1 
                ind_wrong = setdiff(1:self.Npos, [self.reconstruct_ind{:}]);
                probe{i}(:,:,ind_wrong) = repmat(mean(probe{i},3),1,1,length(ind_wrong));
            end
         end
    end
    if par.probe_modes > par.Nrec
        %  orthogonalization of incoherent probe modes 
        if is_used(par, 'fly_scan')
            probe_tmp = probe;
            % orthogonalize the modes with all the other shifted modes 
            for i = 1:par.Nrec
                dx = median(modes{i}.probe_positions - modes{1}.probe_positions);
                probe_tmp{i} = imshift_fft(probe_tmp{i}, dx); 
            end
            probe_tmp = ortho_modes(probe_tmp);  % perform othogonalization 
            probe(1+par.Nrec:par.probe_modes) = probe_tmp(1+par.Nrec:par.probe_modes);
        else
            ind = [1,1+par.Nrec:par.probe_modes];  % skip polyvave/multilayer probe_tmp 
            probe(ind) = ortho_modes_eig(probe(ind));  %% slightly better  
        end
    end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% PREPARE OBJECT, MULTILAYER OBJECT %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    %% updated illumination
    aprobe2 = abs(self.probe{1}(:,:,1)).^2; 
    for ll = 1:par.Nscans
        if par.share_object
            ind = [self.reconstruct_ind{:}];
        else
            ind = self.reconstruct_ind{ll};
        end
        [cache.oROI_s{1}] = find_reconstruction_ROI( modes{1}.probe_positions,self.Np_o, self.Np_p);
        % avoid oscilations by adding momentum term 
        illum_sum_0{ll} = Ggather(set_views(Gzeros(self.Np_o), Garray(aprobe2), 1,1, ind, cache));
    end
    
    
    %% multilayer extension 
    % if object has only a single layer, expand it to multiple using
    % unwrapping 
    if size(self.object,2) > par.Nlayers
        for ll = 1:par.Nscans
            self.object{ll,1} = prod(cat(3,self.object{ll,:}),3); 
        end
        self.object(:,2:end) = []; 
    end
    
    if size(self.object,2) < par.Nlayers
        N_add = par.Nlayers - size(self.object,2); 
        for ll = 1:size(self.object,1)
            obj{ll} = self.object(ll,:); 
            for ii = 1:N_add
                if mod(ii, 2) == 1
                    obj{ll}{end+1} = ones(self.Np_o, 'single') + 1e-9i*randn(self.Np_o, 'single'); % add empty slice at the end 
                else
                    obj{ll}(2:end+1) = obj{ll}; 
                    obj{ll}{1} = ones(self.Np_o, 'single') + 1e-9i*randn(self.Np_o, 'single'); % add empty slice at the beginning 
                end
            end
        end
        self.object = cat(1, obj{:});
    end
    
    % if object has more layers but only one is needed 
    if size(self.object,2) > 1 && par.Nlayers == 1
       for ll = 1:par.Nscans
           object{ll,1} = prod(cat(3,self.object{ll,:}),3); 
       end
       self.object = object; 
    end

    
    for j = 1:par.Nlayers  % add extra layers 
        for i = 1:max(1, par.Nscans * ~par.share_object)
            try
                object{i,j} = self.object{min(end,i),j};
                object{i,j}(1);
            catch
                %% add fully transparent slice at the end 
                object{i,j} = ones(self.Np_o, 'single');
                if size(self.object,2) == 1
                    % swap order of the new layers to keep the original
                    % reconstruction in center 
                    object(i,:) = object(i,end:-1:1); 
                end
            end
            
        end
    end
    
    
    
    for i = 1:numel(object)
        object{i} = single(object{i});
        object{i} =  complex(object{i});
    end
    for i = 1:numel(probe)
        probe{i} = single(probe{i});
        probe{i}  = complex(probe{i});
    end
   
    
    %% STORE RESULTS TO SELF CLASS 
    self.object = object; 
    self.probe = probe; 
    self.modes = modes; 
    self.diffraction = Diffraction; 
    self.noise = Noise; 
    self.mask = Mask; 
    self.background = reshape(Background,1,1,[]);        

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% PRECALCULATE USEFUL VALUES %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
    if ~isfield(self, 'probe_evolution' ) 
        % initial coefficients for OPRP approximation 
        self.probe_evolution(:,1) = ones(self.Npos,1);  % first mode is constant 
    end
    new_probe_modes_ind = 1+(size(self.probe_evolution,2):par.variable_probe_modes); 
    self.probe_evolution(:,new_probe_modes_ind) = 1e-6*randn(self.Npos,length(new_probe_modes_ind));
           
    if par.variable_probe 
        pnorm = norm2(self.probe{1});
        self.probe{1}(:,:,:,2:end) = self.probe{1}(:,:,:,2:end) ./ pnorm(1,1,:,2:end);
        self.probe_evolution(:,2:end) = self.probe_evolution(:,2:end) .* squeeze(mean(pnorm(1,1,:,2:end),3))';
    end
    

    if par.background_detection || ~isempty(self.background)
        %% auto-estimate background correction distribution 
        if isempty(self.mask)
            mask = 0;
        else
            mask = self.mask; 
        end
        if par.background_detection
            background_weight = sum(self.diffraction.^2,3) ./ max(1,sum(~mask,3)); 
            background_weight = imgaussfilt(background_weight,1);
            background_weight = 1./sqrt(max(1e-3, background_weight)) .* ~any(mask,3);
            background_weight = max(0,background_weight - 0.3*mean(background_weight(:)));
            cache.background_weight = ( fftshift_2D(background_weight / sum2(background_weight)));
        end
        
        if isinf(par.background_width)
            cache.background_profile = 1;
        else 
            mdiffr = Garray(fftshift(mean(get_modulus(self,cache,1:self.Npos,false).^2,3)));

            W = par.background_width;
            X = (-self.Np_p(1):self.Np_p(1)-1);
            Y = (-self.Np_p(2):self.Np_p(2)-1);
            [X,Y] = meshgrid(X,Y);

            background_profile = exp(-sqrt( (X/W(1)).^2 +(Y/W(1)).^2)); 
            background_profile = conv2(mdiffr,background_profile, 'same');
            background_profile = background_profile / max2(background_profile);
            
            background_profile = utils.crop_pad(background_profile,self.Np_p); 
            cache.background_profile = gather(fftshift(background_profile)); 
 
        end
         
        if ~isempty(self.diffraction_deform_matrix)
            apply_deform = @(x,D)single(reshape(full(D * double(reshape(x,[],size(x,3)))), size(x)));
            % apply deformation effects caused by tilted sample , be sure to enforce the mask before
            % interpolation, the hotpixels can spread around after the correction
            if isscalar(cache.background_profile)
                cache.background_profile = ones(self.Np_p, 'single'); 
            end
            cache.background_profile = apply_deform(cache.background_profile, self.diffraction_deform_matrix');
        end
    else
    	% looks like cache.background_profile_weight = .. was a typo;
         cache.background_weight = 1; 
    end
    


    for  ll = 1:par.Nscans
        illum_sum_0{ll} = Ggather(illum_sum_0{ll}); 
        cache.MAX_ILLUM(ll) = max(illum_sum_0{ll}(:));
        cache.illum_sum_0{ll} = illum_sum_0{ll}; 
    end

    %% precalculate illumination ROIs
    cache = precalculate_ROI(self,cache, Ggather(sqrt(aprobe2)));

    %% prepare mask needed for subpixel shifts of object views 
    cache.apodwin = single(0.1+0.9*tukeywin(self.Np_p(1),0.05) .* tukeywin(self.Np_p(2), 0.05)');

    if par.initial_probe_rescaling
        %% initial rescaling of probe intensity , just a very rough guess     
        if isfield(par, 'rmvac') && par.rmvac % vac layer removed then modes{end} not a fft, by Zhen Chen
            temp = fft2_safe(self.probe{1}(:,:,1));
            mean_aPsi = mean2(abs(temp.^2)); % bug fixed by Zhen Chen, previous no ^2
        else
            mean_aPsi = mean2(abs(fwd_fourier_proj(self.probe{1}(:,:,1), self.modes{end})).^2); 
        end
        mean_diffraction_intensity = mean(mean2(self.diffraction(:,:,randi(self.Npos, [10,1])).^2));  % take roughly average intensity  % bug fixed by Zhen Chen, previous no ^2

        for ii = 1:par.probe_modes
            self.probe{ii} = self.probe{ii} * sqrt( mean_diffraction_intensity / mean_aPsi); 
        end   
    end

end
