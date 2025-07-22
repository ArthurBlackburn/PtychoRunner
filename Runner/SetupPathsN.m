function exec_params = SetupPathsN(force_local_repo)
% Setup paths for running reconconstructions:
%
% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2021
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
% 
% *******************************************************************

repo_name = 'ps-shelves-test01';
ptychoshelves_subdir = 'PtychoShelvesMS';

path_entry = @(repo_path, varargin) [fullfile(repo_path, repo_name, varargin{:}),';'];

repo_paths =  @(repo_path) [...
             path_entry(repo_path, ptychoshelves_subdir ,'ptycho','utils_EM'),...
             path_entry(repo_path, ptychoshelves_subdir ,'ptycho','utils'),...
             path_entry(repo_path, ptychoshelves_subdir ,'ptycho'),...
             path_entry(repo_path, ptychoshelves_subdir)];

if nargin < 1             
    force_local_repo = false;
end

% for now just a simple way to distinguish between windows and linux.
archstr = computer;
if contains(lower(archstr),'win')
    [~, aa] = system('set COMPUTERNAME');
    if contains(aa,'EXCALIBUR')
            one_loc = 'U:\';
            use_pi_box = false;
            use_lsmql = true;
            show_pibox_recon_window = 0;
            copy_to_monitor_dir = true;
            show_last_recon = false;
            recipe_exec_char = 'a'; 
            plot_results_every = 0;
            repo_path = 'D:\Data\Repos'; 
            if force_local_repo
                addpath(repo_paths(repo_path),'-begin');
            end
            datasets_path = 'C:\Users\Arthur\OneDrive - University of Victoria\EMG\Ptycho\DemoAndTestSets\DatasetFiles';
        elseif contains(aa,'BILBO')
            one_loc = 'C:\Users\stell_000\OneDrive\PtychoShare\';
            use_pi_box = false;
            use_lsmql = true;
            show_pibox_recon_window = 0;
            copy_to_monitor_dir = true;
            show_last_recon = false;
            recipe_exec_char = 'a';
            plot_results_every = 0;
            repo_path = 'C:\Data\Repos';
            if force_local_repo
                addpath(repo_paths(repo_path),'-begin');
            end
            datasets_path = 'G:\Shared Projects\Ptycho\ExampleDataSets\records';
        elseif contains(aa,'HOBO1')
            one_loc = 'S:\';
            use_pi_box = true;
            use_lsmql = false;
            show_pibox_recon_window = 1;
            copy_to_monitor_dir = true;
            show_last_recon = false;
            recipe_exec_char = 'a';
            plot_results_every = 0;
            repo_path = 'U:';
            if force_local_repo
                addpath(repo_paths(repo_path),'-begin');
            end   
        elseif contains(aa,'VAGABOND')    
            one_loc = 'U:\';
            use_pi_box = false;
            use_lsmql = true;
            show_pibox_recon_window = 0;
            copy_to_monitor_dir = true;
            show_last_recon = false;
            recipe_exec_char = 'a';
            plot_results_every = 0;
            repo_path = 'L:';
            if force_local_repo
                addpath(repo_paths(repo_path),'-begin');
            end
            datasets_path = 'G:\Shared Projects\Ptycho\ExampleDataSets\records';
        elseif contains(aa,'WAYFARER')    
            one_loc = 'U:\';
            use_pi_box = false;
            use_lsmql = true;
            show_pibox_recon_window = 0;
            copy_to_monitor_dir = true;
            show_last_recon = false;
            recipe_exec_char = 'a';
            plot_results_every = 0;
            repo_path = 'L:';
            if force_local_repo
                addpath(repo_paths(repo_path),'-begin');
            end
            datasets_path = 'G:\Shared Projects\Ptycho\ExampleDataSets\records';
        else
    end
    base_dir_fieldname = 'base_dir';
        
    exec_params.ptycho_matlab_path = fullfile(repo_path, repo_name, ptychoshelves_subdir,'ptycho');
    exec_params.cSAXS_matlab_path  = fullfile(repo_path, repo_name, ptychoshelves_subdir);
    
%     
else
    
    [~,sysname] = system('hostname');
    if any(contains(sysname,{'eroica','emperor','pastorale'}))
        % Then we're running on Arbutus, but let's point to main storage 
        % location on cedar:
        project_path = '/home/ablackbu/project-spaces/cedar';
    else
        % Elsewhere on main DRA computers point directly into project
        % space:
        project_path = '/project/6042988';
        
    end
    
    % where the experimental data, reconstruction databases, and codebase
    % respectively are stored:
    share_path = fullfile(project_path, 'share');
    pty_path = fullfile(project_path, 'pty');
    repo_path = fullfile(share_path,'repos');
    
    path_extra = repo_paths(repo_path);
    addpath(path_extra);
    
    datasets_path = fullfile(pty_path,'ReconProjects','SummaryDBsLowKVPaper');
    one_loc = fullfile(pty_path,'ReconProjects');
    
    use_pi_box = false;
    use_lsmql = true;
    show_pibox_recon_window = 0;
    copy_to_monitor_dir = false;
    show_last_recon = false;
    recipe_exec_char = 'd';
    base_dir_fieldname = 'base_dir';

    exec_params.ptycho_matlab_path = fullfile(repo_path, repo_name, ptychoshelves_subdir,'ptycho');
    exec_params.cSAXS_matlab_path  = fullfile(repo_path, repo_name, ptychoshelves_subdir);

    plot_results_every = 0; % no plotting
    
    exec_params.cc.project_path = project_path;
    exec_params.cc.pty_path = pty_path;
    exec_params.cc.repo_path = repo_path;
    exec_params.cc.path_extra = path_extra;
    exec_params.cc.datasets_path = datasets_path;
    
end

% if executing as batch file you can also specify:
% -- batch_expt_nums 
% -- log_file

exec_params.one_loc = one_loc; 
exec_params.use_pi_box = use_pi_box;
exec_params.use_lsmql = use_lsmql;
exec_params.show_pibox_recon_window = show_pibox_recon_window;
exec_params.copy_to_monitor_dir = copy_to_monitor_dir;
exec_params.show_last_recon = show_last_recon;
exec_params.recipe_exec_char = recipe_exec_char;
exec_params.plot_results_every = plot_results_every;
exec_params.repo_path = repo_path;
exec_params.base_dir_fieldname = base_dir_fieldname;
exec_params.parallel = false;
exec_params.num_recipes_per_GPU = 1;
exec_params.datasets_path = datasets_path;

exec_params.set_time = now;
