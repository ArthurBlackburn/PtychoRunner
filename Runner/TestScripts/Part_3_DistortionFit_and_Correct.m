% Provides an example of how to determine and correct for pincusion distortion from reconstruction of 
% randomnly oriented gold particles. This was as applied to the manuscript. Steps are provided
% below.
% 
% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2024
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
% 
% *******************************************************************
%  
% ------  Step 1 -----
% Load data from the reconstruction, e.g. results from raw data reconstruction:
% Example file names for the first final reconstruction, using the raw diffraction data without distortion correction:
%     source_reconstruction_file = '...\ReconStore\AuaC_uncorr_R_01_recipe006.mat';
%
% Example file name for a reconstruction performed on the 1-time corrected diffraction data:
%     source_reconstruction_file = '...\ReconStore\AuaC_cubic_corr1_setA_R_01_recipe004.mat';
% OR we could use a later reconstruction:
%     source_reconstruction_file = '...\ReconStore\AuaC_cubic_corr1_R_01_recipe006.mat';
%
% Here we can use the reconstruction from one of the subsets, as only really small differences appear
% in the peak positions in the FT in the last few hundred iterations, so we can look at an earlier 
% recipe:
%
%
% At this point you can make some quick tests to demonstrate fits:
%   source_reconstruction_file = '...\ReconStore\AuaC_cubic_corr1_R_01_recipe006.mat';
%   test1= load(source_reconstruction_file);
%
%   imFT_pattern = GetDistortionFromDiffFit(test1.recon_dat, 'do_plots', [ 0 1 0 1 0 0 0],'fit_method','pattern')
% 
% Which should produce as output (give or take a bit due to random seeding of pattern search) >>
% struct with fields:
% 
%                A: 0.949898910522461
%                k: 0.058350799037398
%              map: [1×1 struct]
%    sf_normalizer: 2.001629307162587
%       A_inv_simp: 1.051146897457186
%       k_inv_simp: 0.056919608007748
%
% imFT_ga = GetDistortionFromDiffFit(test1.recon_dat,  'do_plots', [ 0 1 0 1 0 0 0],'fit_method','ga')
%
% Which should produce as output (give or take a bit due to random seeding of ga search) >>
% struct with fields:
% 
%                A: 0.948424648400910
%                k: 0.068306589932046
%              map: [1×1 struct]
%    sf_normalizer: 2.001629307162587
%       A_inv_simp: 1.052270064018278
%      k_inv_simp: 0.065556730533086
% 
% Let's go forward with the g.a. method:

% First we need to make sure k is determined in the same way as it is applied in lensdistort, i.e. 
% with normalization to maxR.  

centre_index = @(nr)((nr-rem(nr,2))/2)+1;
dpSize = test1.recon_dat.recon_param.asize;
recon_pix_size = test1.recon_dat.recon_pixel_size;
M = dpSize(1); N = dpSize(2); % in all practical instances these are the same...
centre = centre_index(dpSize(end:-1:1)); % order has been swapped to 'x,y' rather 'r,c' as required by lensdistort. 


maxR = sqrt((max(N-centre(1), centre(1)-1))^2 + (max(M-centre(2), centre(2)-1))^2);
reciprocal_pixel_size = 1/(recon_pix_size*N);
max_SF = reciprocal_pixel_size*maxR/1e10; % Put into per Angstrom.

[imFT, fullft] = GetDistortionFromDiffFit(test1.recon_dat,  ...
    'do_plots', [ 0 1 0 1 1 1 0],...
    'fit_method','ga',...
    'sf_normalizer',max_SF);
% you might want to compare a different method, e.g. patternsearch:
% [imFTpa, fullftpa] = GetDistortionFromDiffFit(test1.recon_dat,  ...
%     'do_plots', [ 0 1 0 1 0 0 0],...
%     'fit_method','pattern',...
%     'sf_normalizer',max_SF);
%
% Note that fullft is used in the next block, to illustrate effect on FT of recon.

% Lets work out our reconstructed pixel size....
lens_options = {'borderType','fix',...
                'interpMethod', 'cubic',... % THE INERPOLATION METHOD AFFECTS RECON. QUALITY *
                'padMethod','fill',...
                'padValue',0,...
                'fType',4,...
                'intensCorr',true,...
                'centre',centre}; % Here we use the actual centre of the FT, not the centre that lensdistort uses.
[~, scaling] = lensdistort(ones(512), imFT.k, lens_options{:});
% We could also work out new pixel size without calling the above, as long we are
% absolutely sure to use 'borderType','fix'which sets scale = 1, otherwise the scaling applied
% by lensdistort also has an effect which needs to be considered. (Scaling assumed to be square in
% these cases)
% Also note calibration includes the effect of A i.e. we don't rescale the image but we do correct
% he pxiel size. If A_inv < 1 our DPs angular range was initially too big, or recon pix size was too
% small
new_recon_pix_size = test1.recon_dat.recon_pixel_size * (1/imFT.A_inv_simp) * (1/scaling(1));

fprintf(['Old Recon Pix Size %3.3f pm\n',...
         'New Recon Pix Size %3.3f pm\n'],...
         test1.recon_dat.recon_pixel_size/1e-12, new_recon_pix_size/1e-12);
     
% for first recon: 
% Old Recon Pix Size 23.583 pm
% New Recon Pix Size 24.936 pm     

% for second recon:
% Old Recon Pix Size 24.940 pm (% I rounded to in Excel file)
% New Recon Pix Size 23.806 pm

%% STEP 2 - Inspect the effect on some example diffraction data.
% However, don't get too hung up if things look a bit odd at this stage in FT. This is partly to be
% expected as the reconstruction was performed on distorted diffraction data...
% Do we want to show some test images?
show_test_images = true;

if show_test_images
% lets try this on a test image. This also lets us confirm scaling applied...
    test_image = ones(dpSize);

    lens_options = {'borderType','fix',...
                    'interpMethod', 'cubic',...
                    'padMethod','fill',...
                    'padValue',0,...
                    'fType',4,...
                    'intensCorr',true,...
                    'centre',[]}; % For these 'test' images, it's okay just go the centre of the field of view.
    
    undistorted_image = lensdistort(test_image, imFT.k, lens_options{:});

    figure;

    niceimagesc = @(x) imagesc(x, math.sp_quantile(x,[1e-2, 1-1e-2],10)); 

    subplot(1,3,1);
    imagesc(undistorted_image);
    axis image;
    title('Distorted Ones Image');
    
    counts_in_source = sum(test_image,'all');
    counts_in_distorted = sum(undistorted_image,'all');
    fprintf(['Test image:\nCounts in source image: %3.2f\n',...
             'Count in (un)distorted image %3.2f\nRatio: %3.2f\n'],...
             counts_in_source, counts_in_distorted, counts_in_source/counts_in_distorted);
    
    subplot(1,3,2);
    niceimagesc(abs(fullft.image)); 
    axis image;
    title('Uncorrected FT of Reconstruction');


    subplot(1,3,3);
    corrected_ft_image = lensdistort(abs(fullft.image), imFT.k, lens_options{:});
    niceimagesc(abs(corrected_ft_image)); 
    axis image;
    title('Corrected FT of Reconstruction');

end

%% STEP 3
% If everything looks reasonable, proceed to undistorting the diffraction data to make a new set
% that is passed to tbe next reconstruction:

lens_options = {'borderType','fix',...
                'interpMethod', 'cubic',... % THE INERPOLATION METHOD AFFECTS RECON. QUALITY *
                'padMethod','fill',...
                'padValue',0,...
                'fType',4,...
                'intensCorr',true,...
                'centre',centre}; % Here we use the actual centre of the FT, not the centre that lensdistort uses.
% *Cubic appears to give best results -- from limited experiments -- but has pecularity that the
% interpolation can dip a bit below zero... I think this is dealt with in PS by clipping at zero.

% The location of the stack/s we want to correct:
storage_path = 'D:\Data\Arthur\20250624_CC_Copy1';
save_path = 'D:\Data\Arthur\20250624_testcorrections';
correction_number = 2;

% The file names:
fnames = {'Au_on_aC_full_uncorr.mat',...
          'Au_on_aC_subset_A_uncorr.mat',...
          'Au_on_aC_subset_B_uncorr.mat'};

switch correction_number
    case 1
        corrected_suffix = 'cubic_corr1';
        del_last_chars = 6; % when making new filenames, delete the last 6 chcaracters, and append the suffix above.
        pincushion_coef = imFT.k; % If this first correction on raw data, we can just use forward coefficient directly.
%         sess_savename = ['fitting_session_',datestr(now,'yyyy_mm_dd_HH_MM_SS'),'.mat'];
        sess_savename = ['DistortionFit_',corrected_suffix,'.mat'];
    case 2
        step1_corr = load(fullfile(save_path,'DistortionFit_cubic_corr1.mat'));
        corrected_suffix = 'cubic_corr2';
        del_last_chars = 6; % when making new filenames, delete the last 6 chcaracters, and append the suffix above.
        % For the second correction, pincusion on pincushion, 4th order
        % term can be neglected and the two quadratic terms add. We do this
        % as we don't want computional distort source data twice and have
        % double interpolation errors.
        pincushion_coef = imFT.k + step1_corr.imFT.k; 
%         sess_savename = ['fitting_session_',datestr(now,'yyyy_mm_dd_HH_MM_SS'),'.mat'];
        sess_savename = ['DistortionFit_',corrected_suffix,'.mat'];
    otherwise
        error('Not yet defined...');
end

% Before get carried away with the below let's save the key results of our fitting session which will
% allow us to do the below again without doing another fit etc:

save(fullfile(save_path, sess_savename),...
    'source_reconstruction_file','imFT','recon_pix_size','new_recon_pix_size','dpSize');

%% STEP 4 - If you really want to do it, run the block below! It takes considerable time, a good
% few minutes, to undistort a ~10,000 512 by 512 frames, then save them in compressed v7.3 (hdf5)
% format.
% apply_corrections = true;
apply_corrections = false;

if apply_corrections
% Now lets correct our actual data:
    for fi = 1:length(fnames)
        lname = fullfile(storage_path, fnames{fi});
        load(lname);
        for ii = 1:size(dp,3)
            dp(:,:,ii) = lensdistort(squeeze(dp(:,:,ii)), pincushion_coef, lens_options{:});
        end
        [~,noextname,~] = fileparts(fnames{fi});
        savename = fullfile(save_path,[noextname(1:(end-del_last_chars)), corrected_suffix, '.mat']);
        save(savename,'dp','-v7.3');    
    end
end