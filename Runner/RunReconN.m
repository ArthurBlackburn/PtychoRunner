function RunReconN(runspec)
% Sets up and runs a reconstruction that uses Ptychoshelves, using
% parameters and conditions stored within an spreadsheet (such as created
% within Excel, LibreOffice, OpenOffice, Google Sheets, etc.).
%
% RunReconF(runspec)
% 
% runspec is a struct with fields:
%
% 'datasets_path': the path where information on the datasets, i.e. the
% spreadsheet with 'ExcelFileName', is stored. Some run information is also
% placed here during execution.
%
% 'ExcelFileName': The name of the spreadsheet (*.xls) containing the run
% information.
%
% 'ExcelShortname': A shortname for the spreadsheet that might be appended
% to some produced data output files.
%
% 'ExperimentSheet': The name of the spreadsheet tab that contains
%  information on the Datasets that are to be processed.
%
% 'RecipeSheet': The name of the spreadsheet tab that contains
%  information on the reconstruction recipes that are to be applied to the
%  datasets.
%
% 'RunOnSamples': If the ExperimentSheet contains information on multiple
% datasets that are named in the 'sample_type' column, specifying that name
% here means reconstructions will only be performed on these named datasets.
%
% 'exec_char': Recipes which are to be run are specified with either a single
% char or a string of chars. The letter corresponds to the
% characters placed in the 'exec' column of the RecipeSheet. If there is a
% match the recipe will be executed.
%
% 
% 'FileSaveOverride': A strcuture that allows the portions of the output
% file names specfied in the ExcelSheet to be overwritten with information
% specified here: This struct in turn has fields. It operates by searching
% for a section of the excel specifed output filename and replacing it with
% something specfied here. The fields of the struct are:
%    'search': The segment to search for.
%    'replace_fmt': The format of what to replace the segment with, e.g. if
%    you want to just replace it with a simple string this would just be
%    '%s;
%    'replace_num': Is fed into the above format string. After formatting 
%    this will replace the segment (if the searched-for segment
%    is present), e.g. 
%             struct('search','Rx',...
%                    'replace_fmt','RR_%02.0f',...
%                    'replace_num',2)));
%
% 
% Examples:
% 
% Simple locally stored case:
%
%     RunReconF(struct(...
%     'ExcelFileName','LocallyStoredExample.xls',...
%     'ExcelShortName','LocalTest',...
%     'ExperimentSheet','Datasets',...
%     'RecipeSheet','Recipes',...
%     'RunOnSamples','Au on C',...
%     'exec_char','a'));
%
% Remote execution:
%
%
%     RunReconF(struct(...
%     'datasets_path','/project/user_number/Projects/ReconDBs',...
%     'ExcelFileName','TestRun.xls',...
%     'ExcelShortName','TestR',...
%     'ExperimentSheet','Datasets',...
%     'RecipeSheet','Recipes',...
%     'RunOnSamples','Au on C',...
%     'exec_char','abc',...
%     'FileSaveOverride',struct('search','Rx',...
%                               'replace_fmt','RR_%02.0f',...
%                               'replace_num',2)));
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

% Setup default paths as appropriate for the machine.
exec_params = SetupPathsN(true);

% Default locations for information on datasets:
if isfield(runspec,'datasets_path')
    % possible override:
    datasets_path = runspec.datasets_path;
else
    datasets_path = exec_params.datasets_path;
end


% WHICH EXPERIMENTS / RECIPES
FileLocation = datasets_path;
ExcelFileName = runspec.ExcelFileName;
ExcelShortName = runspec.ExcelShortName;
ExperimentSheet = runspec.ExperimentSheet;
RecipeSheet =  runspec.RecipeSheet;

RunOnSamples = runspec.RunOnSamples;

if isfield(runspec,'exec_char')
    %p possible override:
    exec_params.recipe_exec_char = runspec.exec_char;
end


[experi_def, set_params, recipe_variations] = ...
    LoadFromSpreadSheet(FileLocation, ExcelFileName, ExperimentSheet, RecipeSheet, ...
    RunOnSamples);

experi_def.runspec = runspec;

% When each recipe is executed, the full result is always saved. Below is
% an option to save an additional summmary file containing the key results of
% the all the recipes executed in one batch. This can save tracking down
% files for later analysis if for example the batch was a sweep of
% parameters etc.
savesummary_info = true;

% Name of the summary file name:
SaveNameFunc = @(experinum, recipe_numS) GenSName(ExcelShortName, FileLocation,ExperimentSheet,experinum,RecipeSheet,recipe_numS);

% How to set up the name of the field (or variable name) for a particular run in the summary file:
VarNameFunc = @(experinum, recipe_num) sprintf('%s_%d_%s_%d',ExperimentSheet,experinum,RecipeSheet,recipe_num);

% Whether to run multiple recipes at once using multiple processes. If
% there's only one GPU carefully watch performance and memory requirements
% to see if it is worth it or not.

if isfield(exec_params,'cc') && (length(runspec.exec_char)>1)
    exec_params.parallel = true;
    num_workers = length(runspec.exec_char);
    labletters = runspec.exec_char;
    gpuIndices = 1:length(runspec.exec_char);
else
    exec_params.parallel = false;
end    

% Note our local testing machines only have one GPU. Currently there
% appears to be no advantage to trying to running reconstructions in
% parallel on the same physical GPU; the GPU is usually 'maxed out' with
% just one reconstruction.

%     num_worker = gpuDeviceCount;
%     num_workers = exec_params.num_recipes_per_GPU * gpuDeviceCount;

% When running in parallel, my batch method looks for letters assigned to
% different experiments to decide what runs where. Recipes can have have
% multiple assignments e.g. a, b, c, d, but recipes currently one 1. There
% is lots of room for improvement here!

% if multiple GPUs, how to allocated them to each process:
%     gpuIndices = repmat(1:gpuDeviceCount,exec_params.num_recipes_per_GPU,1);    
%     gpuIndices = gpuIndices(:)';
%     gpuIndices = [1 2 3 4];
% gpuIndices = [1 2];
% gpuIndices = [1 1];


exec_params.verbose_level = 4; % 4 = full verbose level. For decription of other levels see ptychoshelves docs.


% Do parallel execution if specified.
if exec_params.parallel
    poolobj = gcp('nocreate'); % Just look at the existing pool (if there is one). 
    % A pool might have been created before this function has been
    % called. If so we can save time by not recreating it if it's the
    % rights size. We also don't want to create another pool.
    if isempty(poolobj)
        parpool('local',num_workers);
    else
        if poolobj.NumWorkers ~= num_workers
            delete(poolobj);
            parpool('local',num_workers);            
        end
    end
    % Do parallel execution using spmd:
    spmd
        exec_params_local = exec_params;
        exec_params_local.recipe_exec_char = labletters(labindex);
        gpuDevice(gpuIndices(labindex));
        try
            CallerOuter(experi_def, set_params, recipe_variations, exec_params_local, savesummary_info, SaveNameFunc, VarNameFunc);
        catch me
            warning('Error: %s',getReport(me));
        end
    end
else
    CallerOuter(experi_def, set_params, recipe_variations, exec_params, savesummary_info, SaveNameFunc, VarNameFunc);
end

end

function savefilename = GenSName(ExcelShortName , FileLocation,ExperimentSheet,experinum,RecipeSheet,recipe_numS)
    % The name of the summary file, the mat file which contains summary info of all recons done in this batch. 
    if length(recipe_numS) > 1
        savefilename = sprintf('%s_%s_%d_%s_%d_to_%d',ExcelShortName, ExperimentSheet,experinum,RecipeSheet,recipe_numS(1),recipe_numS(end));
    else
        savefilename = sprintf('%s_%s_%d_%s_%d',ExcelShortName, ExperimentSheet,experinum,RecipeSheet,recipe_numS(1));
    end
    savefilename = fullfile(FileLocation , savefilename);
end

function CallerOuter(experi_def, set_params, recipe_variations, exec_params, savesummary_info, SaveNameFunc, VarNameFunc)
    [experi_def_N, recipe_variations_N, experi_num_inds] =  get_recipes(experi_def, set_params, recipe_variations, exec_params);
    experinum = experi_num_inds(1); % For now just doing one experiment run at at time, but would be relatively trivial to loop the body below if needed.
    recipe_numS = experi_def_N.recipe_table.use_variations;
    savefilename = SaveNameFunc(experinum, recipe_numS);
    
    if isfield(experi_def.runspec,'FileSaveOverride')
        searchstring = experi_def.runspec.FileSaveOverride.search;
        repstring = sprintf(experi_def.runspec.FileSaveOverride.replace_fmt, experi_def.runspec.FileSaveOverride.replace_num);
        for exnum =1:find(~ismissing(set_params.save_name),1,'last')
            if ~ismissing(set_params.save_name(exnum)) && ~ismissing(set_params.recon_savepath(exnum))
                ini_filename = set_params.save_name{exnum};
                ini_dirname = set_params.recon_savepath{exnum};
                set_params.save_name{exnum} = replace(ini_filename, searchstring, repstring);
                set_params.recon_savepath{exnum} = replace(ini_dirname, searchstring, repstring);
                for fnamechar = 'ABCD'
                    set_params = customize_filename(experi_def, set_params, fnamechar, exnum);
                end
            end
        end
%         ini_filename = set_params.save_name{experinum};
%         ini_dirname = set_params.recon_savepath{experinum};
%         set_params.save_name{experinum} = replace(ini_filename, searchstring, repstring);
%         set_params.recon_savepath{experinum} = replace(ini_dirname, searchstring, repstring);
    end
    
    % Does this filepath exist? If not, sort it owt gov'nor!
    if ~exist(set_params.recon_savepath{experinum},'dir')
        mkdir(set_params.recon_savepath{experinum});
    end
    
    
    % -- could make this parfor in single exec...
    for recipe_num = recipe_numS
        try
          CallerInner(experi_def_N, recipe_variations_N, set_params, exec_params,...
              recipe_num, experinum, savesummary_info, savefilename, VarNameFunc)
         catch me
            warning('Error: %s',getReport(me));
        end
    end
end

function CallerInner(experi_def_N, recipe_variations_N, set_params, exec_params, ...
    recipe_num, experinum, savesummary_info, savefilename, VarNameFunc)

t = getCurrentTask();
if ~isempty(t)
    labnumber = t.ID;
else
    labnumber = [];
end

set_params2 = ProbeAndObjNames(recipe_variations_N, set_params, recipe_num, experinum, 1);
fprintf('\n\n**** RECIPE NUM %d *** \n\n',recipe_num);
experi_def_N.recipe_table.use_variations = recipe_num; % just to one recipe.
ps_int = PtychoProcessFromTableF_MS(exec_params, experi_def_N, set_params2, experinum, recipe_variations_N);

if savesummary_info
    recon_dat = ps_int.recon_dat;
    if isfield(recon_dat.shelves_op,'cumulative_phase')
        recon_dat.cumulative_phase = ps_int.recon_dat.shelves_op.cumulative_phase;
    end
    if isfield(recon_dat, 'shelves_op')
        recon_dat = rmfield(recon_dat,'shelves_op'); 
        % as this field duplicates a lot of stuff held elsewhere in struct
    end
    varname = VarNameFunc(experinum, recipe_num);
    eval([varname,' = recon_dat']);
    if ~isempty(labnumber)
        sname = sprintf('%s_L%d.mat',savefilename,labnumber);
    else
        sname = [savefilename,'.mat'];
    end    

    runC = ~exist(sname,'file');
    if runC; final_arg = '-v7.3'; else; final_arg = '-append'; end
    % runC = runC+1; 
    % in a simple non parallel execution mode, runC could be set to zero 
    % outside the loop, and condition above would be runC == 0, but as here
    % operating in parallel save a filename with the taskID in it to
    % make it unique.

    save(sname, varname, final_arg);
    clear(varname);
end
            
end

function set_params = customize_filename(experi_def, set_params, fnamechar, exnum)
    rspec_fname = ['replace_filename',fnamechar];
    search_fname = ['CUSTOM_FILENAME_',fnamechar];
    if isfield(experi_def.runspec,rspec_fname) && ~isempty(experi_def.runspec.(rspec_fname))
        set_params.coords{exnum} = replace(set_params.coords{exnum},...
            search_fname, experi_def.runspec.(rspec_fname));
        set_params.file_name{exnum} = replace(set_params.file_name{exnum},...
            search_fname, experi_def.runspec.(rspec_fname));
    end
end

