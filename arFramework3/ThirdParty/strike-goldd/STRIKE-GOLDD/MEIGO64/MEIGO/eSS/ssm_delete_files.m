% $Header: svn://172.19.32.13/trunk/AMIGO2R2016/Kernel/OPT_solvers/eSS/ssm_delete_files.m 770 2013-08-06 09:41:45Z attila $
function delete_files(solver,c_U)
switch solver
    case 'nomad'
        delete fobj_nomad.m
        delete fobj_nomad_omega.m
        delete fobj_nomad_x0.m

        if ~isempty(c_U)
            delete fobj_nomad_Param.m
        end

    case 'n2fb'
        delete objf_n2fb.m

    case 'dn2fb'
        delete objf_dn2fb.m

    case 'fsqp'
        delete fobj_fsqp.m
        if ~isempty(c_U)
            delete constr_fsqp.m
        end

    case 'ipopt'
        delete ipopt_f.m

        if ~isempty(c_U)
            delete ipopt_c.m
        end
end