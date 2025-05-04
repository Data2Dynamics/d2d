% These data structs actually contain no fitted data
nuisance_sets = [arFindData('steadystate_steadystate_phh'), arFindData('steadystate_steadystate_dcf'), arFindData('steadystate_steadystate_dcf_w_smad'), arFindData('qpcr_qpcr_all', 'exact'), arFindData('wb_data_STAT3_TC', 'exact')];

for a = 1 : numel(ar.p)
    [used, m] = arFindParameterUse(a, 0);
    if ( ar.qFit(a) == 1 )
        if used == 0
            ar.p(a) = 0;
            ar.qFit(a) = 2;
            ar.pLabel{a}
        end
        if ( isempty( setdiff( m.di, nuisance_sets ) ) )
            ar.p(a) = 0;
            ar.qFit(a) = 2;
            ar.pLabel{a}
        end
    end
end

% Fix one parameter as scale
sID = arFindPar('scale_socs3_qpcr_qpcr_qpcr_all_nExpID13');
ar.qFit(sID) = 2;
ar.p(sID) = 0;