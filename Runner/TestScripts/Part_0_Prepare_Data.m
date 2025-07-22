% This script demonstrates the pre-processing of the 'raw' data obtained in HDF5 format from Hitachi 
% High-Tech's Azorus software into a Matlab v7.3 *.mat file (also HDF5 based) format that is used 
% in the subsequent files (Part_1and2..., Part_3... etc. This script takes the raw data collected
% for the paper given in the attribution below, recentres it, then forms subsets from it for 
% recnnstructions as described below for the gold on amorphous carbon dataset. 
% 
% Reconstructions are then performed on this data (Part_1and2...), and distortion fitting and 
% correction is carried out (Part3_...). Further reconstructions are then performed on the distortion
% corrected data (Part_1and2). The reconstruction and distortion fitting process is then repeated 
% as described in the paper. Following this the pincushion distortion coefficient is used to correct 
% for distortion in the Au on MoS2 data, as also shown below in this file (Part_0...). Analysis of 
% the results of the reconstructions on the experimental and synthetic data is carried out in Part_4
% and Part_5.
% 
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

data_store_location = 'D:\YourProjects\NC_Dataset\Source\Experimental';
source_file = 'Au_on_aC_source';
full_source_filename = fullfile(data_store_location, [source_file,'.hp']);
%%
% This block prepares a set that is 90 by 90 points, as oppossed to the orignal 114 by 114 mesh. 
% 90 by 90 fits comfortable on the GPU and computer were using at this point.
%
% PrepareData also recentres the data, and puts it in the right order to make it convenient for
% ptychoshelves. The centre position was picked by hand here. It is also possible to recentre in
% Azorus so that it does not have to be done here, in which case the code here would be much, much
% faster.

PrepareData('input_file_or_dir',full_source_filename,...
'positions_trans',struct('trans_coords',5,'rot_angle_deg',0),...
'centre_pos',[203, 214],... % this centre vector was originally chosen manually using graphical selection provided within this function.
'pattern_size_out',[90 90],...
'determine_masks',false,...
'outputDPSizes',{512});

% This produces 4 new files: 
% 
% 1.  Au_on_aC_114by114_90by90_512_dp.mat
% 2.  Au_on_aC_114by114_90by90_512_position.hdf5
% 3.  Au_on_aC_114by114_COORDS.mat
% 4.  Au_on_aC_114by114_rc_512.mat

% 1 - Contains the diffraction data of the 90 by 90 subset.
% 2 - The positional data, in hdf5 format. To get the actual distance over the image plane in nm
%     mulitply these coordinates by 0.1 * scale_in_ps column in the excel file 'Datasets' sheet. For 
%     this work this is 0.1*7.93 = 0.793.
% 3 - The 'voltage coordinates' and some related parameters extracted directly from the *.hp file.
%     Note that these voltage coordinates, will require some scaling transposing, rotating, etc to
%     get to real space. The field positions_trans.trans_coords give the transforamtion pair used.
% 4 - Contains all the recentred diffraction data. Note here as the input size is 512 by 512 and the
%     output size in 512 by 512, then inevitably some collected data will be cropped out of the field.

% For the purpose of placing the files the repository, the files were renamed:
% 1 -> Au_on_aC_full_uncorr.mat
% 2 -> Au_on_aC_full_posn.hdf5
% 3 and 4 were not placed in the repository, as they are not directly used to produced
% reconstructions presented in the paper. _However_, #3 and #4 are used as intermediates in the
% production of sub-sets A and B below:  

% Now we make subsets A and B from the full set.
% The subsets are taken from a checkerboard pattern of the orignal grid.
% Note that if the original grid had an uneven dimensions, then some rows or columns of the 
% resulting subset mesh will be +/- 1 in length than adjacent rows / or columns. Not a problem for
% ptycho, just a 'beware' in case math somewhere depends on number in rows or columns.

% Prepare the checkerboard mask:
full_coords = load(fullfile(data_store_location,[source_file,'_COORDS.mat']));
scan_size = full_coords.pattern_shape;
mask_cols = false(scan_size);
mask_cols(:,2:2:end) = true;
mask_rows = false(scan_size);
mask_rows(2:2:end,:) = true;
mask = xor(mask_cols, mask_rows);
% mask was arranged so that cols represent x and row reps y...
% scan proceeds in a raster type fashion, so transpose to get sequence:
mask = mask'; % transpose.
data_mask_linear = mask(:);

%%
newfull = load(fullfile(data_store_location,[source_file,'_rc_512.mat']));

% Make set A diffraction data:
dp = newfull.dp(:,:, data_mask_linear);
save('set_A_new.mat','dp','-v7.3');

% Make set A positions:
setposns = full_coords.voltage_coords(data_mask_linear,:);
output_name = fullfile(data_store_location,'set_A_new_posns.hdf5');
transpose_type = full_coords.positions_trans.trans_coords;
rot_angle_deg = full_coords.positions_trans.rot_angle_deg;
% Due to  earlier code a factor 10 is introduced for Angstrom defaults..
convert_to_angstrom = true;
nm_V = 1;

gen_coords_for_ptychoshelves(setposns, output_name, nm_V, rot_angle_deg, convert_to_angstrom, transpose_type);
%%
% Make set B diffraction data. Note the tilde ~ which negates the masks:
dp = newfull.dp(:,:, ~data_mask_linear);
save('set_B_new.mat','dp','-v7.3');

% Make set A positions:
setposns = full_coords.voltage_coords( ~data_mask_linear,:);
output_name = fullfile(data_store_location,'set_B_new_posns.hdf5');

gen_coords_for_ptychoshelves(setposns, output_name, nm_V, rot_angle_deg, convert_to_angstrom, transpose_type);

%%
% As for the MoS2 datasets
source_file = 'MoS2_source';
full_source_filename = fullfile(data_store_location, [source_file,'.hp']);

PrepareData('input_file_or_dir',full_source_filename,...
'positions_trans',struct('trans_coords',5,'rot_angle_deg',0),...
'centre_pos',[232, 231],... % this centre vector was originally chosen manually using graphical selection provided within this function.
'pattern_size_out',[45 45],...
'determine_masks',false,...
'outputDPSizes',{512});

% take distortion coef from Au on aC (prepared in part 3):
step1_corr = load(fullfile(data_store_location,'DistortionFit_cubic_corr1.mat'));
step2_corr = load(fullfile(data_store_location,'DistortionFit_cubic_corr2.mat'));
% seems to be scaled approx by camera length and cube of mag... values used first time round are 
% a little different to those here, probably due to g.a. solver giving a slightly different value
% due to random seed, or other small offsets in code. (23.8 and 16.57 are pixel sizes, see xls).
pincushion_coef = (step2_corr.imFT.k + step1_corr.imFT.k)*((23.8/(16.57*2))^3);
% Experimentally it appears that scaling as ^3 gives a reasonable fit, as seen by straight line MoS2
% pattern. ^3 has some plausibility due to displacements from Cs varying with magnification^3, but
% really need to do more work (another paper) to see if this holds all the time and is supported by 
% theory. Likely depends on details of projector system... For now though we get a more or less 
% straight line on MoS2 data diffraction data (i.e. no visible pincushion), so reconstruction 
% quality should be vastlk improved, and certainly plausible. Main test and demo of a resolution 
% measure comes from FRC and Au on C data, which has been more carefully and accurately calibrated in
% comparison to Au on MoS2 data. 

apply_corrections = true;

centre_index = @(nr)((nr-rem(nr,2))/2)+1;
dpSize = [512 512];
centre = centre_index(dpSize(end:-1:1)); % order has been swapped to 'x,y' rather 'r,c' as required by lensdistort. 


lens_options = {'borderType','fix',...
                'interpMethod', 'cubic',... % THE INERPOLATION METHOD AFFECTS RECON. QUALITY *
                'padMethod','fill',...
                'padValue',0,...
                'fType',4,...
                'intensCorr',true,...
                'centre',centre}; % Here we use the actual centre of the FT, not the centre that lensdistort uses.

if apply_corrections
% Now lets correct our actual data:
    lname = fullfile(data_store_location, 'MoS2_source_45by45_512_dp.mat');
    load(lname);
    for ii = 1:size(dp,3)
        dp(:,:,ii) = lensdistort(squeeze(dp(:,:,ii)), pincushion_coef, lens_options{:});
    end
    [~,noextname,~] = fileparts(lname);
    savename = fullfile(data_store_location,[noextname,'_corr.mat']);
    save(savename,'dp','-v7.3');    
end
%%
% now pad out to 1024 by 1024:
dp_src = dp;
dp = utils.crop_pad(dp,[1024, 1024,size(dp_src,3)]);
savename = fullfile(data_store_location,[noextname,'_1024_corr.mat']);
save(savename,'dp','-v7.3');  



