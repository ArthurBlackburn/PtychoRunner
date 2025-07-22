function [experi_def_N, recipe_variations_N, run_inds_N] =  get_recipes(experi_def, set_params, recipe_variations, exec_params)
% Used to extract which recipe numbers to run. Used by PtychoProcessFromTableF_MS.m
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

    %Define which experiments and recipes to run:
    experi_def_N = experi_def;
    recipe_variations_N = recipe_variations;
    if experi_def.summary_table.use_ys_to_define_experis
        run_sets = strcmpi(set_params.do_recon, exec_params.recipe_exec_char);
        sample_sets = contains(set_params.sample_type,experi_def.summary_table.run_on_samples);
        run_inds_N = find(run_sets & sample_sets)';
    end

    % Define which recipes to run:
    if ~isfield(experi_def.recipe_table,'use_variations') || isempty(experi_def.recipe_table.use_variations)
        experi_def_N.recipe_table.use_variations = find(contains(recipe_variations.exec,exec_params.recipe_exec_char))';
    end
    % and sort out the recipe sub-steps
    if iscell(recipe_variations.sub_steps)
        for ii = 1:length(recipe_variations.sub_steps)
            recipe_variations_N.sub_steps{ii} = str2num(recipe_variations.sub_steps{ii}); %#ok<ST2NM>
        end
    else
        recipe_variations_N.sub_steps = num2cell(recipe_variations.sub_steps);
    end


end