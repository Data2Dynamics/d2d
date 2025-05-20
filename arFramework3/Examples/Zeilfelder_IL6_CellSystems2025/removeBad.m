try
    % There were a few experimental measurements of SOCS3 mRNA which were
    % clear outliers prior to model analysis. These are removed.
    bad = arFindData('nExpID24', 'input_dcf', 500, 'input_il6', 10, 'input_noggin', 0, 'state', 'SOCS3RNA');
    badTimes = find(ar.model.data(bad).yExp(:, strcmp(ar.model.data(bad).y, 'SOCS3_qpcr_scaled'))<-2);
    ar.model.data(bad).yExp(badTimes, :) = NaN;
    if numel( badTimes ) > 0
        disp('Removed bad SOCS3 mRNA outlier');
    end
catch
end

try
    % There were a few experimental measurements of SMAD mRNA which were
    % clear outliers prior to model analysis. These are removed.
    bad = arFindData('nExpID24', 'input_dcf', 500, 'input_noggin', 0, 'state', 'SMAD6RNA');
    badTimes = find(ar.model.data(bad).yExp(:, strcmp(ar.model.data(bad).y, 'SMAD6_qpcr_scaled'))<-1.1);
    ar.model.data(bad).yExp(badTimes, :) = NaN;
    if numel( badTimes ) > 0
        disp('Removed bad SMAD6 mRNA outlier');
    end
catch
    
end

try
    % This experiment contains no BMP response and is removed.
    badsmad6 = arFindData('qpcr_qpcr_APAP_nExpID5');
    
    for a = 1 : numel(badsmad6)
        ar.model.data(badsmad6(a)).qFit(strcmp(ar.model.data(badsmad6(a)).y, 'SMAD6_qpcr_scaled'))=0;
    end
    
    ar.qFit(arFindPar('scale_smad6_qpcr_qpcr_qpcr_APAP_nExpID5')) = 0;
    
    disp( 'Removed bad SMAD6 mRNA time course' );
catch
    
end

fprintf('Don''t forget to run arLink!\n');