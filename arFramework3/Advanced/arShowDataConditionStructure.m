function arShowDataConditionStructure

global ar

for jm=1:length(ar.model)
    fprintf('\nmodel #%i: %s\n',jm,ar.model(jm).name);
    for jd=1:length(ar.model(jm).data)
        fprintf('\n\tdata #%i: %s',jd,ar.model(jm).data(jd).name);
        
        jc = ar.model(jm).data(jd).cLink;
        fprintf(' -> condition #%i',jc);
        nss = sum(ar.model(jm).condition(jc).qSteadyState);
        if(nss>0)
            fprintf(' [steady-states=%i]',nss);
        end
        fprintf('\n');
        
        for jc=1:length(ar.model(jm).data(jd).condition)
            fprintf('\t%s = %s\t',ar.model(jm).data(jd).condition(jc).parameter, ...
                ar.model(jm).data(jd).condition(jc).value);
        end
        fprintf('\n');
    end
end
    

