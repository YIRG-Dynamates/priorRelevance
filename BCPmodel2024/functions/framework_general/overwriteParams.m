function param_settings = overwriteParams(params,param_settings,fit_settings)
%Overwrite the param_settings with the 'params' values for fitted params

if ~isempty(params)
    for i=1:fit_settings.num_conds
        for j=1:numel(fit_settings.fit_param_nrs_per_cond{i})
            param_nr = fit_settings.fit_param_nrs_per_cond{i}(j);
            param_name = fit_settings.fit_param_names{param_nr};
            param_settings.(param_name)(i) = params(param_nr);
        end
    end
end

end %[EoF]
