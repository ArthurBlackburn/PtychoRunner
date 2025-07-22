function [imFT, ft, plot_data] = ...
    GetDistortionFromDiffFit(recon_dat,  varargin)
% Takes an image / reconstruction of gold particles, takes the FFT,
% extracts a radial profile, and fits of this to the modelled Au diffraction pattern. 
% The model assumes that the collected diffraction pattern is a pincushion distorted version of
% the original.
%
% Parameters are then extracted from this that can be used to correct the FFT
% and - importantly - the source diffraction patterns that were used to
% create the reconstruction in the first place. 
%
% Example - results from raw data reconstruction:
% test1= load('D:\Data\ReconJune2025\AuaC_AsColl_R_01_recipe006.mat');
% imFT_pattern = GetDistortionFromDiffFit(test1.recon_dat,  'do_plots', [ 0 1 0 1 0 0 0],'fit_method','pattern')
%
% imFT_ga = GetDistortionFromDiffFit(test1.recon_dat,  'do_plots', [ 0 1 0 1 0 0 0],'fit_method','ga')
% imFT_ga = GetDistortionFromDiffFit(test1.recon_dat,  'do_plots', [ 1 1 1 1 1 1 1],'fit_method','ga')
%
% Pattern search and ga give slightly different values. Trying to understand which is better fit at
% this point is of limited use: instead lets do another reconstruction with some values of
% distortion correction applied, that makes this at least a bit better! Then we can do another
% reconstruction, which will have sharper peaks (as the diffraction data is then more physically
% consistent with a physical sample). Once we have sharper peaks, then it's perhaps more use
% wondering which is the best method.
%
% e.g.
% test2= load('D:\Data\ReconJune2025\AuaC_first_corr_R_01_recipe006.mat');
% im2FT_pattern = GetDistortionFromDiffFit(test2.recon_dat,  'do_plots', [ 0 1 0 1 0 0 0],'fit_method','ga')
%
% Looks like g.a. lines up peaks better at higher SF imo.
%
% General idea:
% load up reconstructed image, get the diffraction pattern, in radial form,
% then find fitting parameters that give best correlation with the
% known diffraction pattern for gold.
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

do_tests = true;

% Just an example of 'reconstruction data', which must contain the (complex) object moodel:
if isempty(recon_dat)
    source_recon_name = 'D:\Data\ReconJune2025\AuaC_AsColl_R_01_recipe006.mat';
    load(source_recon_name,'recon_dat');  
end

p = inputParser;
% G: The maximum spatial frequency in the data. Best left at default (i.e. set by assiging []), 
% unless the max freq in your experimental data is way beyond where there is actually any useful
% information.
p.addParameter('sf_normalizer',[]);
% The model data file for performing the fit upon:
p.addParameter('model_file','AuModelDiffDataExt_20kV.mat'); % This default should be located in your path.
% fp: An index into the diffraction model diffraction which is after the direct beam (at 0), but before
% first diffraction peak, so that normalization can be performed on the strongest diffraction peak 
% rather than the direct beam.
p.addParameter('fp',50); 
% Type of window to apply before taking the FT:
p.addParameter('window_type','hamming', @(x)ismember(lower(x), {'hamming', 'hann'}));
p.addParameter('custom_edgesize',[]);
p.addParameter('fit_method','anneal',@(x) ismember(x,{'ga','anneal','fmin','pattern','brute'}) );
p.addParameter('corr_metric', 1, @(x) isnumeric(x) );
p.addParameter('do_plots',true,  @(x) isnumeric(x) || islogical(x))
if isstruct(recon_dat)
    p.addParameter('lambda',recon_dat.lambda);
    p.addParameter('recon_pix_size',recon_dat.recon_pixel_size);    
elseif ismatrix(recon_dat)
    p.addParameter('lambda',1); % over-ride above or set if not present
    p.addParameter('recon_pix_size',1);    
end

p.parse(varargin{:});
params = p.Results;

do_plots = params.do_plots;
if isscalar(do_plots)
    do_plots = params.do_plots & true(1,6);
end
% note on plots:
% 1: The (2d) FT of the experimental data.
% 2: The radial average of the experimental and the model data.
% 3: Just a test to illustrate an example mapping
% 4: Experimental data and distorted perfect model
% 5: Normalized Forward and Backward Mappings
% 6: Comparison of full inverse and simplified inverse mapping
% 7: Illustration of Mappings


% STEP 1: Get the FT of the reconstruction:
[res, ft ] = radial_FT_from_recon(recon_dat, 'ave_meth','sums');
res.radial_ave_111pk_normed = res.radial_ave/(max(res.radial_ave(params.fp:end)));

% loads model, where col 1 is 1/d (A) and col 2 is expect intens.
% See L:\em_utils\Ptycho\Processing\DistortFit.m on how to create this
% data.
modl = load(params.model_file);
fname = fieldnames(modl);
modl = modl.(fname{1});   


%%
if do_plots(2)
    figure;
    plot_data.experi_radial.x = res.radial_ord(params.fp:end);
    plot_data.experi_radial.y = res.radial_ave_111pk_normed(params.fp:end);
    plot(plot_data.experi_radial.x, plot_data.experi_radial.y);    
    xlabel('d-spacing, 1/A');
    ylabel('Peak Normalized Intensity');
    legcell = cell(1);
    legcell(1) = {'Reconstruction Profile'};
    hold on;
    plot_data.model_radial.x = modl(:,1);
    plot_data.model_radial.y = modl(:,2);
    plot(plot_data.model_radial.x, plot_data.model_radial.y,'r');
    legcell(end+1) = {'Model Profile'};
    legend(legcell);
    title('Diffraction Radial Spectra: Reconstruction and Model, Before Fitting');
end

%%
% Now lets' find a fit of reconstruction to model, using pincushion distortion model. Here
% undistorted radial data (r) is mapped to distorted space radial coordinate, s:
%
% s = A*r.*(1+k.*(r.^2));
%
% A is the scale factor;
% K is the pincushion coefficient (see also DistortFit.m)

% Below is an example of what this looks like with 'typical' coefficient values:

if do_plots(3) && do_tests
    A = 0.97;
    % Note: k must be positive. Putting a -ve here is not quite the same as barrel, see also DistortFit.m for inverse!
    k = 0.56; 
    s_fwd = @(r, A, k) A*r.*(1+k.*(r.^2));
    map_func = @(r) s_fwd(r, A, k);
    title_string = 'Forward Mapping Test';
    figure;
    axm = axes;
    plot_mapping(map_func, title_string, axm)
end


if isempty(params.sf_normalizer)
    params.sf_normalizer = res.radial_ord(end); % The first column contains spatial freqs. End element is max spatial freqeuncy.
end

% Prepare the fit function
fit_func = @(pvec) DistortFit(pvec(1), pvec(2), ...           % A, k
    res.radial_ord(params.fp:end),              ...    % Experimental x
    res.radial_ave_111pk_normed(params.fp:end), ...    % Experimental y
    params.sf_normalizer, ...     % The spatial frequency in x to use as normalizer  
    params.model_file, params.corr_metric);
lb = [0.5 0]; ub = [1.5 0.7];

% ---------- DO THE FITTING! --------
% Here I'm simply minimizing an error metric ('difference signal') between the model and
% reconstruction. Other more advanced methods, would involve peak fitting and identification, which
% might be more reliable... but also would require some user interaction and likely a significantly
% greater programming effort. So simply minimizing an error signal is used here for now, and seems 
% to work well provided you keep checking the results.


switch params.fit_method
    case 'ga'
        options = optimoptions('ga');
        options.FunctionTolerance = 1e-7;
        options.MaxGenerations = 500;
        best_p = ga(fit_func,2,[],[],[],[],lb, ub,[],options); 
    case 'anneal'
        rng default % For reproducibility, reset Matlab's random number generator here.
        k_guess = 0.08*params.sf_normalizer^2;    % from some experiments a reasonable starting guess:
        A_guess = 0.95;
        p_guess = [A_guess, k_guess];
        best_p = simulannealbnd(fit_func,p_guess,lb,ub);
    case 'fmin'
        k_guess = 0.08*params.sf_normalizer ^2;
        A_guess = 0.95;
        p_guess = [A_guess, k_guess];
        options = optimoptions('fminunc','Display','off',...
            'Algorithm','quasi-newton');
        best_p = fminunc(fit_func,p_guess,options);
    case 'pattern'
        options = optimoptions('patternsearch');
        options.CompletePoll = 'on';
        options.CompleteSearch = 'on';
        k_guess = 0.08*params.sf_normalizer^2;
        A_guess = 0.95;
        p_guess = [A_guess, k_guess];
        best_p = patternsearch(fit_func,p_guess,[],[],[],[],lb,ub,[],options);
    case 'brute'
        % Calculate error for a full grid of A and k, and take the minimum.
        %
        % For circularized 2D data, 'brute force' and using error_metric = ...
        % -sum(model_diff_intens_distorted .* experi_data_y .* model_d_pts_norm.^2);
        % seems to work best.!
        a_range = 0.9:0.005:1.1;
        k_range = 0:0.005:0.4;
        [a_mesh, k_mesh] = meshgrid(a_range, k_range);
        p_vals = [a_mesh(:), k_mesh(:)];
        dfs = zeros(size(a_mesh));
        for ii = 1:size(p_vals,1)
            dfs(ii) = fit_func(p_vals(ii,:));
        end
        [~, best_ind] = min(dfs(:));
        best_p = p_vals(best_ind,:);
    otherwise
        error('Unknown model fitting specifer')
end

% ---------------------------
% The resulting 'best fit' vector contains A and k.
imFT.A = best_p(1);
imFT.k = best_p(2);

% So now we have an A and k that takes model to the experimental data.
% The map and sf_normalizer are:
[~, imFT.map, imFT.sf_normalizer] = fit_func(best_p);

% These are the points from the experiment that were mapped into the
% model (ideal) space:
experi_pts = res.radial_ord(params.fp:end);
experi_pts_norm = experi_pts/imFT.sf_normalizer;
% So in the model (ideal space) the radial coords were found from the inverse map
model_pts_norm = imFT.map.inv(experi_pts_norm, imFT.A, imFT.k);

% And so the points below from the model (found at model_pts_norm) whne plotted at experi points
% show the distorted model and experimental data:
model_diff_intens_distorted = interp1(modl(:,1)/imFT.sf_normalizer, modl(:,2), model_pts_norm,'spline');

if do_plots(4)
    figure;
    plot(experi_pts, model_diff_intens_distorted,'m'); hold on;
    plot(experi_pts, res.radial_ave_111pk_normed(params.fp:end),'k');
    title('Experimental data and distorted perfect model (after fitting)');
end
% If things have gone well, the 2 plots should have similarly positioned peaks. If not tweek fit
% params and try again...

% Lets now look at the forward and reverse mappings:
r_d = 0.01:0.01:1;
mapped_fwd = imFT.map.fwd(r_d, imFT.A, imFT.k);
mapped_backward = imFT.map.inv(r_d, imFT.A, imFT.k);

if do_plots(5)
    figure;
    plot(r_d, mapped_fwd,'b');
    hold on;
    plot(r_d, mapped_backward,'r');
    legend({'Forward','Backward'});
    title('Normalized Forward and Backward Mappings');
end

% The accurate inverse mapping is in fact a complicated function (see DistortFit.m), but we can do
% a further fit on this to provide a simplified mapping, using something simpler:
% Spefically:
%       s = r.*(1./(1+k.*(r.^2))); % Note the 1./ over...

imFT.map.r_inv_simp = @(r, A_inv, k_inv) A_inv*r.*(1./(1+k_inv.*(r.^2)));

% Set up fittype and options.
[xData, yData] = prepareCurveData( r_d, mapped_backward ); % This the accurate model going distorted space to undistorted space.

ftm = fittype( 'A*x.*(1./(1+k.*(x.^2)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.5 0.5];
% Fit model to data.
[fitresult, ~] = fit( xData, yData, ftm, opts );


imFT.A_inv_simp = fitresult.A;
imFT.k_inv_simp = fitresult.k;

if do_plots(6)
    figure;
    plot(r_d, mapped_backward,'x');
    hold on;
    plot(r_d,fitresult(r_d),'r');
    legend({'Full model','Fitted curve'});
    title('Comparison of full inverse and simplified inverse mapping')
end

% So in most instances can use this form to map to undistorted form and save a bit of compuational time.
% However, for most image undistort algorithms it is just the forward mapping that we need.

% Nonetheless the above is handy as it can allow us to work out the new reconstructed pixel size,
% quite easily, remembering that in the case of Au the initial calibration was done of the 111 ring, 
% so we can map that back to undistorted space and find out how many pixels that takes up now ->
% then scale accordingly.

if do_plots(7)
    figure;
    ax1 = subplot(1,3,1);
    fwd_map = @(r) imFT.map.fwd(r, imFT.A, imFT.k);
    plot_mapping(fwd_map, 'Forward Map', ax1);
    ax2 = subplot(1,3,2);
    rev_map = @(r) imFT.map.inv(r, imFT.A, imFT.k);
    plot_mapping(rev_map, 'Reverse Map', ax2);
    
    ax3 = subplot(1,3,3);
    rev_map = @(r) imFT.map.r_inv_simp(r, imFT.A_inv_simp, imFT.k_inv_simp);
    plot_mapping(rev_map, 'Simplifed Reverse Map', ax3);
    sgtitle('Illustration of Mappings');
end

end

function plot_mapping(map_func, title_string, ax_hand)
    % Just a quick way to plot out what the mapping looks like.
    x_test = -10:2:10;
    [XX, YY] = meshgrid(x_test);
    [theta, RR] = cart2pol(XX, YY);
    % Normally the below would apply, as in lensdistort
    max_RR = max(RR(:));
    % However, if sf_normalizer or other is taken as the limit of our radial profile which 
    % only extends out to maximum in X, then it would be more representative to use:
    % max_RR = max(XX(:));
    % to give a more accurate visualization of distortion (or at least less reinterpretation by the viewer!)    
    RR = RR/max_RR;
    RR_dist = map_func(RR(:));
    RR_dist = RR_dist*max_RR;
    RR_dist = reshape(RR_dist,size(theta));
    [XX_dist, YY_dist] = pol2cart(theta, RR_dist);
    surf(XX_dist,YY_dist,ones(size(XX_dist)),'FaceColor','None','EdgeColor','r' ); view(2); 
    axis square;
    hold on;
    % ax_hand = surf(XX,YY,ones(size(XX_dist)),'FaceColor','None','EdgeColor','b' );
    set(ax_hand,'GridColor','None');

    quiver(XX, YY, (XX_dist - XX), (YY_dist - YY),1,'AutoScale','off');
    title(title_string);

%     RR_undist = r_inv(RR_dist/max_RR,A,k)*max_RR;
%     [XX_undist, YY_undist] = pol2cart(theta, RR_undist);
%     surf(XX_undist,YY_undist,ones(size(XX_dist)),'FaceColor','None','EdgeColor','g','LineStyle','--');


end