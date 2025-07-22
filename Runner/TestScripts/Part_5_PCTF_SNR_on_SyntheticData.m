% Determine the SNR and PCTF from the ptychographic data reconstructed on 
% synethetic test objects, where the simulated data has the same or very similar
% reconstructed pixel size and similar (or larger) 20 keV electron beam illumination 
% half-angle as in the experimental data.
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

% Where the copy of the publicly available reconstruction data is locally copied:
recon_repository_location = 'D:\YourProjects\NC_Dataset\Reconstructions';

% The filename of the ground truth test image. The ground truth test image was
% designed to have a reasonably broad spectral content and contain  series of 
% double-dot structures that can be used to allow assessment of resolution in the
% classic Rayleigh sense of the meaning of resolution.
ground_truth_location = 'D:\YourProjects\NC_Dataset\Reconstructions';
padded_test_image_filename =   fullfile(recon_repository_location ,'testobject_2048pix_padded.mat');

ground_truth = load(padded_test_image_filename);

fnames = {
    'synth_MoS2conds_20m_R_01_recipe030.mat',... % 20 mrad
    'synth_MoS2conds_10m_R_01_recipe036.mat',... % 10 mrad
    'synth_MoS2conds_5m_R_01_recipe042.mat',... % 5 mrad
    'synth_aCconds_15m_R_01_recipe048.mat' , ... % 15 mrad
    'synth_aCconds_10m_R_01_recipe054.mat' , ... % 10 mrad
    'synth_aCconds_5m_R_01_recipe060.mat' }; ... % 7.5 mrad

titles = {...
'$MoS_{2}: \alpha =$ 20 mrad',...
'$MoS_{2}: \alpha =$ 10 mrad',...
'$MoS_{2}: \alpha =$ 5 mrad',...
'$aC: \alpha =$ 15 mrad',...
'$aC: \alpha =$ 10 mrad',...
'$aC: \alpha =$ 7.5 mrad'};

num_files = numel(fnames);

recon_pix_sizes = [16.57e-12, 16.57e-12, 16.57e-12,...
                   23.80e-12, 23.80e-12, 23.80e-12];
%%
% Regions on the outside of a scan also don't get information coming from
% as many overlapping regions as on the inner region of reconstruction. If
% we trim down the compared region further, by DP_SIZE, then the inner
% region (in terms of over-lapping regions) is approximately translational
% symmetric in integer intervals of scanning lattice basis vectors. So
% trimming down is done to reduce edge effects and give a scan size
% independent measure.
% As DP_SIZE = 512, lets just set compare region to 1024 pixels (1500 - 512, 
% then rounded to 2^(integer) to help speed ffts a bit):
% compare_central_region_size = 1024;

% This is what went into manuscript
compare_central_region_sizes = [2000, 2000, 2000, 1700, 1700, 1700];

% But this would have been more consistent with hindsight...:
% compare_central_region_sizes = 1024*ones(1, num_files);

% Setup for alignment of images:
% Use the phase of objects to align on:
measure = @(x) angle(x);
% extent of the seearch for good alignment
offset_extent = 40;

rc = cell(1, num_files);
%
% Load in, align and trim the data:
for ii = 1:num_files

    % Load the reconstruction:
    tmp = load(fullfile(recon_repository_location , fnames{ii}));
    % account for roation mismatch if any.
    rc{ii}.object = rot90(tmp.recon_dat.obj_int{1},3);

    % Get the ground truth, appropriately trimming the object.
    ground_truth_trimmed = PtyAn.get_trimmed(ground_truth, compare_central_region_sizes(ii));
    rc{ii} = PtyAn.get_trimmed(rc{ii}, compare_central_region_sizes(ii));

    [rc{ii}, ground_truth_trimmed, rc{ii}.best_ind_ABS] = PtyAn.simple_align(rc{ii}, ground_truth_trimmed, [], measure,...
    {'search_step',1,'subsample',2,'offset_extent',offset_extent});    
    rc{ii}.gt = ground_truth_trimmed;
end
%%
% If you want to see the alignments between the reconstruction and GT:
for ii = 1:num_files
    PtyAn.plot_quick(rc{ii}.gt, rc{ii}, 'best_match');
end
%%
%
% Now calculate the PCTF and SNR
for ii = 1:num_files
    % Perform phase leveling. Small misalignments between the images
    % results in a phase ramp, which is of no consequence to SNR and PCTF - it is just an
    % image misregistrations. (This phase ramp is the basis of some
    % sub-pixel alignment routines)
    levelled_fieldname = 'best_match_levelled';
    rc{ii} = PtyAn.phase_level(rc{ii}, rc{ii}.gt.best_match, rc{ii}.best_match, levelled_fieldname);
    
    % Get the PCTF
    % Note that this routine performs phase unwrapping, before making
    % comparison. If there are strong phase changes, etc that do not
    % unwrap well note that the PCTF may be misleading.
    params = {'tftype','pctf', 'rc_name', levelled_fieldname, 'tlevel',0.02,'for_correction',true};
    rc{ii} = PtyAn.get_TFs(rc{ii}.gt.best_match,  rc{ii}, recon_pix_sizes(ii),  params);

    % Now get ready to determine the SNR.
    % As there is no such thing as absolute phase, reconstruction and object may
    % also have a phase piston (shift) between them. Let's remove that before determining SNR:    
    mean_uw_gt_phase = mean(rc{ii}.pctf.unwrapped_gt(:));
    mean_uw_recon_phase = mean(rc{ii}.pctf.unwrapped_obj(:));
    poffset = mean_uw_gt_phase - mean_uw_recon_phase;

    % A bit more subtle is that amplitude can be scaled due to probe
    % updates not always enforcing re-normalization back to the measured beam
    % current. Ptychographic reconstruction forces diffraction patterns to
    % match experiment, but the in general the probe current has soom room to drift: 
    % not many (?) reconstruction algorithms take as input what the measured beam current 
    % / flux actually is, other than (perhaps) taking the counts of the intial probe guess, 
    % and _trying_ to keep the probe intensity at that level.
    % For exmaple, probe intensities (0.5, 1) and object transmissions (2, 1) give exit 
    % intensities (1, 1) so both result could / would 
    % result in the same fourier error metric. However, in the
    % case of amplitude drift, the relative amplitude changes (i.e.
    % contrast) should remain same regardless of amplitude 'drift'...
    %
    % Thus, let's normalized objects before comparison of amplitudes:
    amp_ratio = mean(abs(rc{ii}.gt.best_match)./abs(rc{ii}.best_match),'all');


    % Now let's find the SNR. Note this is a full complex comparison
    rc{ii} = PtyAn.get_SNR_measures(rc{ii}, rc{ii}.gt, recon_pix_sizes(ii), ...
        {'recon_fieldname',levelled_fieldname,...
        'gt_fieldname','best_match',...
        'amp_mult', amp_ratio,...
        'phase_shift',poffset, ...
        'bin_factor',4, ...
        'nbins',50,...
        'compon','complex'});   
end
%%
% Plot the characteristics
% The PCTF:
fhandle0 = [];
for ii = 1:num_files
    fhandle0 = PtyAn.plot_TFs(rc{ii}, {'pctf_simple'}, fhandle0);
end
sgtitle('PCTF');
ylabel('PCTF');
xlabel('Spatial Frequency, m^{-1}');
legend(titles,'Interpreter','latex','Location','SouthWest');
xlim([0.1 3]*1e10); ylim([0 1.8]);
%%
% The SNR
fhandle1 = [];
for ii = 1:num_files
    fhandle1 = PtyAn.plot_TFs(rc{ii}, {'SNR_rad'}, fhandle1, @(y) log10(y));
end

sgtitle('SNR (complex)');
ylabel('Log_{10}(SNR)');
xlabel('Spatial Frequency, m^{-1}');
legend(titles,'Interpreter','latex','Location','NorthEast');

xlim([0.15 2.2]*1e10); ylim([0 5.3]);