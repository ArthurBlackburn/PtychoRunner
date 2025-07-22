function [set_params] = ProbeAndObjNames(recipe_variations_N, set_params, recipe_num, experinum, ...
    name_gen_func)
% Used for setup in PtychoProcessFromTableF_MS.m, RunRecon etc.
% *******************************************************************
% Author & Copyright: Arthur M. Blackburn
% Year              : 2021
% Contact           : ablackbu@uvic.ca / arthur.blackburn@gmail.com
% Citation and Attribution: If possible, please cite the related publication
% given at https://github.com/ArthurBlackburn/PtychoRunner.
% Otherwise, thanks would be appreciated if you find this code useful.
% 
% *******************************************************************

if isempty(name_gen_func)
    % an example function that should be created externally
    engine = 'PSS';
    mat_file_name = 'test';
    name_gen_func = @(experi_num, r_num) ...
        sprintf('%s_%s_r%0.0d_%s',...
         engine, ...
         recipe_variations_N.algorithm{r_num}, ...
         experi_num, ...
         mat_file_name);
end

if name_gen_func == 1
    name_gen_func = @(experi_num, r_num) ...
        sprintf('%s',...
         set_params.save_name{experi_num});    
end

if any(strcmpi(recipe_variations_N.probe_model(recipe_num),{'step_numbered','step_numbered_masked'}))
    if ismember('probe_expt_num',recipe_variations_N.Properties.VariableNames) && ...
            isnumeric(recipe_variations_N.probe_expt_num(recipe_num)) && ...
            recipe_variations_N.probe_expt_num(recipe_num) ~= experinum && ...
            ~recipe_variations_N.force_same_expt_num(recipe_num) 
        
        load_probe_e_num = recipe_variations_N.probe_expt_num(recipe_num);
    else
        load_probe_e_num = experinum;
        recipe_variations_N.probe_expt_num(recipe_num) = experinum;
        % in principle should not need to do the above, probe_load_name
        % should be used everywhere when needing to determine the name of
        % the probe to load, but just in case there are cases whent this is
        % not th case, set it also above to cover the when
        % recipe_variations_N.force_same_expt_num(recipe_num) == true;
    end
    
    set_params.probe_load_name{recipe_num} = ...
        name_gen_func(load_probe_e_num, recipe_variations_N.probe_recipe_num(recipe_num));
    
end

if any(strcmpi(recipe_variations_N.object_model(recipe_num),{'step_numbered','step_numbered_masked'}))
    if ismember('object_expt_num',recipe_variations_N.Properties.VariableNames) && ...
            isnumeric(recipe_variations_N.object_expt_num(recipe_num)) && ...
            recipe_variations_N.object_expt_num(recipe_num) ~= experinum && ...
            ~recipe_variations_N.force_same_expt_num(recipe_num) 
        
        load_object_e_num = recipe_variations_N.object_expt_num(recipe_num);
    else
        load_object_e_num = experinum;
        recipe_variations_N.object_expt_num(recipe_num) = experinum;
    end
    set_params.object_load_name{recipe_num} = ...
        name_gen_func(load_object_e_num, recipe_variations_N.object_recipe_num(recipe_num));
end



end

