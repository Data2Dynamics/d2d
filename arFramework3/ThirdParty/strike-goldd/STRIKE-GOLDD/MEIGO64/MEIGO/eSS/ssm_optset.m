% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ssm_optset.m 770 2013-08-06 09:41:45Z attila $
function [opts] = ssm_optset(default,opts)

%Esto es para poner el campo 'local' el ultimo
if isfield(opts,'local') & not(isempty(opts.local))
    temp=opts.local;
    opts=rmfield(opts,'local');
    opts.local=temp;
else
    opts.local=default.local;
end

if not(isempty(opts));
    opt_names=fieldnames(opts);
    
    opt_names_local=fieldnames(opts.local);
    
    
    default_names=fieldnames(default);
    default_names_local=fieldnames(default.local);
    
    low_opt_names = lower(opt_names);
    low_default_names = lower(default_names);
    
    low_opt_names_local = lower(opt_names_local);
    low_default_names_local = lower(default_names_local);
    
    for i=1:length(opt_names)-1
        j = strmatch(low_opt_names{i,:},low_default_names,'exact');
        if isempty(j)
            error(sprintf(['Unrecognized field name %s'], opt_names{i,:}));
        end
        value = opts.(opt_names{i,:});
        default.(default_names{j,:}) = value;
    end

    for i=1:length(opt_names_local)
        j = strmatch(low_opt_names_local{i,:},low_default_names_local,'exact');
        if isempty(j)
            error(sprintf(['Unrecognized field name %s'], opt_names_local{i,:}));
        end
        value = opts.local.(opt_names_local{i,:});
        default.local.(default_names_local{j,:}) = value;
    end    
    
end
opts=default;
