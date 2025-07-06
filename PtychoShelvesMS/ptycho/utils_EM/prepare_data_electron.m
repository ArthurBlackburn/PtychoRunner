% This script prepares data for PtychoShelves from clean EMPAD dataset,
% *.mat file provided in our paper: 
%  Chen, Z., et al., Electron ptychography achieves atomic-resolution limits set by lattice vibrations. arXiv: 2101.00465.
% Zhen Chen @ Cornell University, 3/27/2021

%% Step 1: download sample data (rawdata.mat) from the link provided in 
% PARADIM website to release, stay tuned: https://data.paradim.org/
% to start a data_dp.mat (3 dimension) is required.

clear;
%% Step 2: load data
data_dir = '<path_for_data_storage>'; %change this
fmat='<filename>.mat';  % *.mat file

scan_number = 01; %Ptychoshelves needs
%%
Np_p = 256; % size of diffraction patterns used during reconstruction. can also pad to 256
ADU=580; % counts per electron 
alpha0=21.4; % mrad
rbf=26; % radius of center disk in pixels
voltage=300; % kev
df=-200; % estimated defocuse
scanstep=0.41; % scan step size in angstrom
rot_ang=0; % relative angle between scan and diffraction
%% Step 3: go back to .../fold_slice/ptycho and pre-process data
Np_p=[Np_p,Np_p];
fdata=fullfile(data_dir,fmat);
load(fdata);
% pad cbed

crop_idx=[1,64,1,64]; % start from smaller data
cbed=cbed(:,:,crop_idx(1):crop_idx(2),crop_idx(3):crop_idx(4));
[ndpy,ndpx,npy,npx]=size(cbed);
if ndpy < Np_p(1) % pad zeros
    dp=padarray(cbed,[(Np_p(1)-ndpy)/2,(Np_p(2)-ndpx)/2,0,0],0,'both');
else
    dp=crop_pad(cbed,Np_p);
end

dp = dp / ADU; % convert to electron count
dp=reshape(dp,Np_p(1),Np_p(2),[]);
Itot=mean(squeeze(sum(sum(dp,1),2))); %need this for normalizting initial probe

% calculate pxiel size (1/A) in diffraction plane
lambda = 12.3986./sqrt((2*511.0+voltage).*voltage);
% [~,lambda]=electronwavelength(voltage);
dk=alpha0/1e3/rbf/lambda; %%% PtychoShelves script needs this %%%

%% Step 4: save CBED in a .hdf5 file (needed by Ptychoshelves)
% scan_number = 1; %Ptychoshelves needs
save_dir = fullfile(data_dir,num2str(scan_number,'%02d'));
mkdir(save_dir)

save(fullfile(save_dir,'data_dp.mat'),'dp','-v7.3');

%% Step 5: prepare initial probe
dx=1/Np_p(1)/dk; %% pixel size in real space (angstrom)
cs = 0;
probe=generateProbeFunction(dx,Np_p(1),0,0,df,cs,1,voltage,alpha0,0);
probe=probe/sqrt(sum(sum(abs(probe.^2))))*sqrt(Itot)/sqrt(Np_p(1)*Np_p(2));
probe=single(probe);
% add parameters for PtychoShelves_electron
p = {};
p.binning = false;
p.detector.binning = false;

save(fullfile(save_dir,'probe_initial.mat'),'probe','p')
copyfile(strcat(mfilename('fullpath'),'.m'),save_dir);

%% prepare probe position 
probe_positions_0=position_generate(npx,npy,scanstep,scanstep, rot_ang);
hdf5write(fullfile(save_dir,'data_position.hdf5'),'/probe_positions_0',probe_positions_0);
