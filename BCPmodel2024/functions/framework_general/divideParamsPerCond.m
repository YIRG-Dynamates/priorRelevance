function param_settings_cond = divideParamsPerCond(param_settings,fit_settings)
% Subdivide the param_settings struct per condition

param_names = fieldnames(param_settings);
num_cond = fit_settings.num_conds;
param_settings_cond = cell(1,num_cond);
for c=1:num_cond
    for j=1:numel(param_names)
        param_settings_cond{c}.(param_names{j}) = param_settings.(param_names{j})(c);
    end
end

end %[EoF]