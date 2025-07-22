% This script re-runs the calculations used to produce the FRC characteristics given in the
% manuscript.
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

leg_titles = {...
    '4: flat', '19: rand', '64: flat + apt', '70: rand + apt';...
    '4: THR', '19: THR', '64: THR', '70: THR'};

color_order = {'r','g','b','m'};
marker_order = {'--','-x','-+','-*'};
figure;

xlabel('Spatial Freq');
ylabel('FRC');


% The parameters below were taken from manual selections from the produced graphs and running 
% particulate_FRC_res manually (i.e. without the offsets provided here, and using 'autoalign'. See 
% particulate_FRC_res). These provide the alignment of reconstructions and ROI of reconstructions. 
% Values that gave a higher resolution limit were taken to indicate better alignment. 
alig_data = ...
{4,  [25 -1], [1081 1024 3233 3316];... % Recipe 04
 19, [4  15], [941  942  3297 3197];... % Recipe 19
 64, [17 29], [1061 1045 3281 3359];... % Recipe 64
 70, [9 -2],  [1026 1034 3229 2963]};   % Recipe 70

recon_store = 'D:\Data\ReconJune2025';
results = cell(1,4);
for ii = 1:4

    
    r0_name = sprintf('AuaC_second_corr_setA_R_01_recipe%03d.mat',alig_data{ii,1});
    r1_name = sprintf('AuaC_second_corr_setB_R_01_recipe%03d.mat',alig_data{ii,1});
    xg1 = alig_data{ii,2}(1);
    yg1 = alig_data{ii,2}(2);
    fov = alig_data{ii,3};


    params ={'passtoalign',...
            {'FOV_corners', fov,...
            'guessx',xg1,'guessy',yg1,...
            'searchForBest',true, ...
            'do_align',false,...
            'plot_level',0,...
            'searchGridPts',8,...
            'autoplotclose', false,...
            'image_prop','complex',...
            'use_window',false,... %% might boost things - and noted that there is no need to taper twice.
            'roisize',[512 512],...
            'window_type','hann',...
            'show_2D_fourier_corr',false}}; 
       
    r0 = load(fullfile(recon_store,r0_name));
    r1 = load(fullfile(recon_store,r1_name));

    results{ii} = particulate_FRC_res(r0, r1, params);

    plot(results{ii}.spatial_freq, results{ii}.FRC, color_order{ii});
    hold on;
    plot(results{ii}.spatial_freq, results{ii}.thresh, [color_order{ii},marker_order{ii}],'MarkerSize',3);
    drawnow;
end

%%
% Get Model and experimental FT
% Model FFT:
model_file = "L:\ps-shelves-test01\Runner\AuModelDiffDataExt_20kV.mat";
AuMod = load(model_file);

% Calculate the experimental FT:
% Load a full-set (as opposed to subset A or B) reconstruction in order to get the radial profile
full_recipe_number = 64;
r_full = load(fullfile(recon_store, sprintf('AuaC_second_corr_R_01_recipe%03d.mat',full_recipe_number)));
res = radial_FT_from_recon(r_full.recon_dat, 'do_plots',true, 'ave_meth', 'sums');
%%
% to replot the above in case figures are closed:
figure;
for ii = 1:4
    plot(results{ii}.spatial_freq, results{ii}.FRC, color_order{ii});
    hold on;
    plot(results{ii}.spatial_freq, results{ii}.thresh, [color_order{ii},marker_order{ii}],'MarkerSize',3);        
end

plot(AuMod.Model(:,1), AuMod.Model(:,2),'b--');
leg_titles{1,end+1} = 'model data';
plot(res.radial_ord(50:end),res.radial_ave(50:end)./max(res.radial_ave(50:end),[],'all'),'r','LineWidth',1.5);
leg_titles{2,end} = 'experimental FT';

legend(leg_titles);
xlabel('Spatial Frequency $\mathrm{(\AA^{-1})}$','Interpreter','latex','Fontname','Arial');
ylabel({'Peak Normalized Intensity or','Fourier Ring Correlation Measure'},'Fontname','Arial');
xlim([0 2.0]);
ylim([0 1.0]);

%%
% Plot for Figure 2 in Manuscript:
% Load the model diffraction data:

F2_titles = {...
    'FRC: Flat Init Phase, No aperture', 'FRC: Random Init Phase & Aperture','FRC: 1/2 Bit Threshold',...
    'Radially Summed FT', 'Model Diffraction Intensity'};

F2_colors = {[1 0 0],[1 0 0],[0 0.5 0],'k',[0 0 1]};
F2_linestyles = {'-','--','-','--','-'};
linewidths = 1.5 * ones(1,5);
F2_new = figure;

plt_nm = 1; 
ii = 1;
plot(results{ii}.spatial_freq, results{ii}.FRC, ...
    'Color', F2_colors{plt_nm}, 'LineStyle', F2_linestyles{plt_nm},...
    'LineWidth',linewidths(plt_nm));
hold on;

plt_nm = 2; 
ii = 4;
plot(results{ii}.spatial_freq, results{ii}.FRC, ...
    'Color', F2_colors{plt_nm}, 'LineStyle', F2_linestyles{plt_nm},...
    'LineWidth',linewidths(plt_nm));

plt_nm = 3;
ii = 4; % Note threshold from all results is the sames, see this by looking at plots above.
plot(results{ii}.spatial_freq, results{ii}.thresh, ...
    'Color', F2_colors{plt_nm}, 'LineStyle', F2_linestyles{plt_nm},...
    'LineWidth',linewidths(plt_nm));

plt_nm = 4;
plot(AuMod.Model(:,1), AuMod.Model(:,2), ...
    'Color', F2_colors{plt_nm}, 'LineStyle', F2_linestyles{plt_nm},...
    'LineWidth',linewidths(plt_nm));

plt_nm = 5;

plot(res.radial_ord(50:end),res.radial_ave(50:end)./max(res.radial_ave(50:end),[],'all'),...
    'Color', F2_colors{plt_nm}, 'LineStyle', F2_linestyles{plt_nm},...
    'LineWidth',linewidths(plt_nm));

legend(F2_titles);
fontname(F2_new,'Arial');
legend boxoff
xlabel('Spatial Frequency $\mathrm{(\AA^{-1})}$','Interpreter','latex','Fontname','Arial');
ylabel({'Peak Normalized Intensity or','Fourier Ring Correlation Measure'},'Fontname','Arial');
xlim([0 2.0]);
ylim([0 1.0]);

%%


