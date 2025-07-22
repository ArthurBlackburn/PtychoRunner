function recipe = check_pty_recipe(recipe)
% A function for setting defaults an doing some checks on recipes for reconstructions.
%
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

default_recipe.params = {...
    'dp_SquareRoot',                { true;  true} ;
    'dp_FFTShift',                  { true;  true};
    'dp_Transpose',                 { false; false};
    'dp_FlipHorizontal',            { false; false};
    'dp_FlipVertical',              { false; false};
    'dp_Pedestals_incoherent',      { 0; 0};
    'dp_Pedestals_coherent',        { 0; 0};
    'probeTransforms_unityOnAdd',   { 0; 0};
    'rp_scaleFactor',               { 1; 1};
    'rp_subPixel',                  { false; false};
    'rp_Object_update',             { false; true};
    'rp_Object_updateFactor',       { 1; 1};
    'rp_Object_phaseOnly',          { false; false};
    'rp_Probe_update',              { true; true};
    'rp_Probe_updateFactor',        { 1; 0.1};
    'rp_Pos_Random_enabled',        { false; false};
    'rp_Pos_acceptValue',           { 1.01 ; 1.01 };
    'rp_Pos_Random_startRange',     { 2;  2};
    'rp_Pos_Random_maxRange',       { 10; 10};
    'rp_Pos_Random_stepFactor',     { 1 ; 1};
    'rp_Pos_Random_stepCount',      { 400; 400};
    'rp_Diff_Drift_useOffset',      { false; false};
    'rp_Diff_Drift_calculateOffset',{ false; false} }';
default_recipe.params = struct(default_recipe.params{:});
default_recipe.step_names = {'Probe normalize'; 'Refinement'};
default_recipe.iterations = [1; 50];
default_recipe.request_results = [true; true];

if ~isstruct(recipe)
    disp('First argument should be struct. As not supplied, using default parameters');
    % create default params
    recipe = default_recipe;
end

if numel(recipe.iterations) < numel(recipe.params)
    fprintf(['Number of iterations for all steps are not specified!\n',...
             'Ignoring steps where iterations were not specified.\n']);
    recipe.params = recipe.params(1:numel(recipe.iterations));
end

% should also throw away steps with zero iterations specified.
if any(recipe.iterations == 0)
    fprintf('Throwing away some steps with zero iterations specified');
    recipe.params = recipe.params(recipe.iterations > 0);
end

%% Check on parameters
default_p_names = fieldnames(default_recipe.params);
supplied_p_names = fieldnames(recipe.params);

% check for existence of parameter fields, and do quick sanity checks on supplied
% data, filling in blanks with defaults when not supplied.
for idx = 1:numel(default_p_names)
    if ~isfield(recipe.params,default_p_names{idx}) 
        for step_ind = 1:numel(recipe.params)               
            if step_ind == 1; use_step = 1; else; use_step = 2; end
            recipe.params(step_ind).(default_p_names{idx}) = default_recipe.params(use_step).(default_p_names{idx});
            fprintf('* Default applied for parameter %s in step %02.0f \n', default_p_names{idx}, step_ind);
        end
    end
end

% check for fields that appear in the input that are do not appear in the
% scheme, and warn against them.
for idx = 1:numel(supplied_p_names)
   if ~any(strcmpi(supplied_p_names{idx},default_p_names))
       fprintf('* Parameter supplied was removed as not in schema: %s\n', supplied_p_names{idx});
       recipe.extra_params.(supplied_p_names{idx}) = recipe.params.(supplied_p_names{idx});
%            recipe.params(step_ind) = rmfield(recipe.params(step_ind), supplied_p_names{idx});
       recipe.params = rmfield(recipe.params, supplied_p_names{idx});
   end
end


%% Check step names
if numel(recipe.step_names) < numel(recipe.params)
    % grow with names 'step_n':
    for grow_ind = (numel(recipe.step_names)+1):numel(recipe.params)
       place_idx = circshift([grow_ind, 1], int8(isrow(recipe.params))); 
       recipe.step_names{place_idx(1), place_idx(2)} = sprintf('STEP_%02.0f',grow_ind);
    end
end

if numel(recipe.request_results) < numel(recipe.params)
    append_elems = numel(recipe.params) - numel(recipe.request_results);
    append_shape = circshift([append_elems, 1], int8(isrow(recipe.params)));
    append_true = ones(append_shape, 'logical');
    if isrow(recipe.params)
        recipe.request_results = [recipe.request_results, append_true];
    else
        recipe.request_results = [recipe.request_results; append_true];
    end
    % Geez there must be a simpler way...
end

end