function [opts] = vns_optset(default,opts)

if not(isempty(opts));
    opt_names=fieldnames(opts);
    default_names=fieldnames(default);
    
    low_opt_names = lower(opt_names);
    low_default_names = lower(default_names);

    for i=1:length(opt_names)
        j = strmatch(low_opt_names{i,:},low_default_names,'exact');
        if isempty(j)
            error(sprintf(['Unrecognized field name %s'], opt_names{i,:}));
        end
        value = opts.(opt_names{i,:});
        default.(default_names{j,:}) = value;
    end
end
opts=default;
