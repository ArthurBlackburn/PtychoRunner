function [corr, map, sf_normalizer] = DistortFit(A, k, experi_data_x, experi_data_y, sf_normalizer, model_file, error_metric)
% Function to provide a quality of fit (correlation - corr, output) between
% experimental polycystalline gold diffraction data and model diffraction data,
% using scaling (A) and quadratic (k) pincushion distortion on the model
% data. Currently set up to use simulated Au 20 kV data by default. Look inside to see detail!
%
% Model data comes from CrystalMaker output, converted for example as below:
%
%{%
% 
% % Default version from Xtalmaker has reflections limit of 0.7 A (1/d = 1.42)
%     model_data_file = 'C:\Users\Arthur\OneDrive\Documents\Papers\JSM2022\Gold Profile 20 kV.txt';
% % Extended version has reflections limit of 0.4 A (1/d = 2.5). Set this in Xtalmakeer options. 
%     model_data_file = 'C:\Users\Arthur\OneDrive\Documents\Papers\JSM2022\Gold Profile Extended 20kV.txt';
% % Read in the model
%     Model = readmatrix(model_data_file);
% % Crystalmakers spits out not enough decimal places in its text output so due to rounding errors, 
% % so there are some useless duplicated data points. Maybe CrystalMaker has some hidden
% % preference setting to increase output decimal places, but I can't find it... 
% % For now thouugh lets remove duplicated points:
%     [~, IA] = unique(Model(:,1));
%     Model = Model(IA,:);
% % Normalize the data to the highest peak:
%     peak_val = max(Model(:,2));
%     Model(:,2) = Model(:,2)/peak_val;
% % See if it makes senses
%     figure; 
%     plot(Model(:,1), Model(:,2));
% % Save it somewhere handy, in mat format so we don't have to do this again:
% % save('L:\em_utils\Basic_Utils\AuModelDiffData_20kV.mat','Model');
%     save('L:\em_utils\Basic_Utils\AuModelDiffDataExt_20kV.mat','Model');
% %
% Or if you have 2D output from Xtalmaker, the function below will ciruclarize it (form a radial
% profile from it). See the related function:
%
%   CircularizeRefList()
%}
%
% For now to see how to use this function, you will have to either look inside, or look at calling
% functions in demos.
% 
% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2022
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
% 
% *******************************************************************


    % The fit method:
    if nargin < 7
        error_metric = 1;
    end
    
    % Model file
    if nargin < 6
        model_file = 'AuModelDiffData_20kV.mat'; % This default should be located in your path.
    end
    
    % The maximum spatial frequency in the data.
    % Best left a default (i.e. set by assiging []), unless the max freq in your experimental data is 
    % way beyond where there is actually any useful information.
    if nargin < 5
        sf_normalizer = [];
    end
    
    % Check on distortion parameters:
    if ~all(size(A) == size(k)) && ~(isscalar(A) || isscalar(k)) % XOR-ish
        error('Size of inputs A and B should either be same size, unless A or k is a single value');
    end
    
    
    num_vals = max(numel(A), numel(k));
    if numel(A) >= numel(k)
        corr = zeros(size(A));
    else
        corr = zeros(size(K));
    end
    
    for ii = 1:num_vals
        
        if isscalar(A)
            Ai = A;
        else
            Ai = A(ii);
        end

        if isscalar(k)
            ki = k;
        else
            ki = k(ii);
        end

        if ii ~= num_vals
            corr(ii) = DistortFit_single(Ai, ki, experi_data_x, experi_data_y, sf_normalizer, model_file, error_metric);
        else
            % On the last run, give the map. In practice we only should request the map if we put
            % one k and on A in. 
            [corr(ii), map, sf_normalizer] = ...
                       DistortFit_single(Ai, ki, experi_data_x, experi_data_y, sf_normalizer, model_file, error_metric);
        end
    end

    % Why am I passing through sf_normalizer?
    % ... Work in progress; TO DO: I am intending to put scaling operations carried out in
    % GetDistortFromDiffFit here.

end

function [err_met, map, sf_normalizer] = DistortFit_single(A, k, experi_data_x, experi_data_y, sf_normalizer, model_file, error_metric)

    persistent model_data
    
    % To save loading data from disc or passing it all the time:
    if ~exist('model_data','var') || isempty(model_data)
        model_data = load(model_file);
        fname = fieldnames(model_data);
        model_data = model_data.(fname{1});    
    end
    
    if nargin < 5 || isempty(sf_normalizer)
        sf_normalizer = model_data(end,1); % The first column contains spatial freqs. End element is max spatial freqeuncy.
    end
    
    % Note here our forward mapping on the diffracted data is
    map.fwd = @(r, A, k) A*r.*(1+k.*(r.^2));
    % However, this forward mapping is not actually used here! I just put it here for reference and
    % completeness. What we actually need is the inverse mapping of the above, which is not so straight
    % forward. Using maple, the inversion (i.e. a rearrangement for r in terms of s, where s = map,fwd) 
    % can be found (My reference working is at 
    % "C:\Users\Arthur\OneDrive\Documents\Papers\JSM2022\distortion_invert.mw")
    %
    % This gives the mapping as:
    map_inv_base = @(s, A, k) -(12.^(1./3)) .* ...
        (-(A .^ (1 ./ 3)) .* (sqrt((27 .* s .^ 2 .* k + 4 .* A .^ 2)) .* sqrt(3) + 9 .* s .* sqrt(k)) .^ (2 ./ 3) ...
        + (A .* 12 .^ (1 ./ 3))) .* ...
        (k .^ (-1 ./ 2)) .* (A .^ (-2 ./ 3)) .* (sqrt((27 .* s .^ 2 .* k + 4 .* A .^ 2)) .* sqrt(3) + 9 .* s .* sqrt(k)) .^ (-1 ./ 3) ./ 6;
    % This might run into issues (0/0 or Infs) if k = 0. To avoid that case causing code to bomb out, or
    % forcing to think about L'hopital's rule etc, let just make K very small number instead of zero.
    % This is used in Lensdistort code and is fine in practiical situations.
    map.inv = @(s, A, k) map_inv_base(s, A, (k==0)*1E-9 + k);

    % Let's use our maximum spatial frequency (sf_normalizer) from the model to normalize for 
    % consistency by defaults. Idea is just to get things in to the range < 1. We should 
    % eventually put a check here that experimental max (1/d) is less than the model max (1/d), 
    % as in this case there would be an unused tail of experimental data that is not fitted to
    % anything.

    model_sf_pts_norm = map.inv(experi_data_x/sf_normalizer, A, k);
    % Here we take our experimental distorted space positions and map that back to r the undistorted
    % space, i.e. r is function of s, so this is the inverse mapping. Why do this: Because our model 
    % data is smooth and noise free (coming from infinite electron dose model) it can interpolated 
    % at any r with with minimal or no errors.
    % In comporisionn trying to interpolate noisy experimental data is likely just going to end up 
    % decreasing our S:N. Better to just map noise free model.r
    model_diff_intens_distorted = ...
        interp1(model_data(:,1)/sf_normalizer, model_data(:,2), model_sf_pts_norm,'spline');

    % Now the tricky bit. What do we use for our error metric. A more thorough investigation will perhaps 
    % one day reveal which one works most reliably. Or maybe some one has already done this??
    % The fit metric below is minimized by the fitting routine:
    
    switch error_metric
        case 1
            % Simple sqrt of the sum of squared errors (SSE):
            err_met = sqrt(sum((model_diff_intens_distorted - experi_data_y).^2,'all'));
        case 2
            % My proposed alternatives see to work better: 
            % SSE maybe puts too much emphasis on lower freq bigger stronger peak?
            % how about linear multiplier to increase effect of upper peaks?:
            % Counter this with a linear ramp (model_d_pts_norm).
            % Note: normalization of SFs is going to have an effect here!
            err_met = sqrt(sum((model_diff_intens_distorted.*model_sf_pts_norm - experi_data_y).^2,'all'));
        case 3
            % or a quadratic ramp:    
            err_met = sqrt(sum((model_diff_intens_distorted.*model_sf_pts_norm.^2 - experi_data_y).^2,'all'));
        case 4
            % product of the two funcs:
            err_met = -sum(model_diff_intens_distorted .* experi_data_y);
        case 5
            % product of the two funcs with linear multiplier.    
            err_met = -sum(model_diff_intens_distorted .* experi_data_y .* model_sf_pts_norm);
        case 6
            % product of the two funcs with quadratic multiplier.
            % (Seems pretty good for sharp peak models such as MoS2)
            err_met = -sum(model_diff_intens_distorted .* experi_data_y .* model_sf_pts_norm.^2);
        otherwise
            error('Unknown fit metric specified');
    end

end