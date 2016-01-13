function PPL_fitScale(m,c,ix,t,scaling_ps,onlyProfile)
global ar
pReset = ar.p;
qFit = ar.qFit;

ps_fix = [];
if(~exist('scaling_ps','var') || isempty(scaling_ps))
    yNames = ar.model(m).data(c).yNames(ix);
    for i=yNames;
       tmp_find = strfind(ar.model(m).data(c).py,yNames{i}); 
       for j=1:length(ar.model(m).data(c).py)
          if(~isempty(tmp_find{j}))
             ps_fix = [ps_fix,  strcmp(ar.pLabel,ar.model(m).data(c).py(j))];
          end
       end
    end
    
elseif(isnumeric(scaling_ps(1)))
    ps_fix = scaling_ps;
elseif(ischar(scaling_ps(1)))
    for i = scaling_ps
        ps_fix = [ps_fix,find(strcmp(scaling_ps(i),ar.pLabel)==1)];
    end
end
 
do_refit_obs = 1;
do_vali = 1;
 
pstmp = [];
if(onlyProfile)
    for it=t
        pstmp = [pstmp; squeeze(ar.model(m).data(c).ppl.ps(it, ix,:,:))];               
        
    end
else
    pstmp = [squeeze(ar.model(m).data(c).ppl.ps_high(:,ix,:)); squeeze(ar.model(m).data(c).ppl.ps_low(:,ix,:))];
end
 
if(ar.model(m).data(c).has_tExp)
    ar.model(m).data(c).yExpUB(:,ix) = ar.model(m).data(c).yExpSimu(:,ix);
    ar.model(m).data(c).yExpLB(:,ix) = ar.model(m).data(c).yExpSimu(:,ix);
else
   ar.model(m).data(c).tExp = ar.model(m).data(c).tFine;
   ar.model(m).data(c).yExp(:,ix) = ar.model(m).data(c).yFineSimu(:,ix);
   ar.model(m).data(c).yExpUB(:,ix) = ar.model(m).data(c).yExpSimu(:,ix);
   ar.model(m).data(c).yExpLB(:,ix) = ar.model(m).data(c).yExpSimu(:,ix);
end

ar.model(m).data(c).yFineUB(:,ix) = ar.model(m).data(c).yFineSimu(:,ix);
ar.model(m).data(c).yFineLB(:,ix) = ar.model(m).data(c).yFineSimu(:,ix);

ar.model(m).data(c).qFit(ix) = 1;
ar.qFit(:) = 2;
ar.qFit(ps_fix) = 1;
 
hbar = waitbar(0, 'Please wait...');
 
for i = 1:size(pstmp,1)
    ar.p(:) = pstmp(i,:);
 
    hbar = waitbar(i/size(pstmp,1), hbar, 'Please wait...');
    try
        arSimu(false, true);
        arSimu(false, false);
        if(do_refit_obs)
            arFit(true);
        end
        arSimu(false, true);
        arSimu(false, false);
 
        
        if(~do_vali)
            if(ar.model(m).data(c).has_tExp)
                q = ar.model(m).data(c).yExpSimu(:,ix) > ar.model(m).data(c).yExpUB(:,ix);
                ar.model(m).data(c).yExpUB(q,ix) = ar.model(m).data(c).yExpSimu(q,ix);
                q = ar.model(m).data(c).yExpSimu(:,ix) < ar.model(m).data(c).yExpLB(:,ix);
                ar.model(m).data(c).yExpLB(q,ix) = ar.model(m).data(c).yExpSimu(q,ix);
            end
            q = ar.model(m).data(c).yFineSimu(:,ix) > ar.model(m).data(c).yFineUB(:,ix);
            ar.model(m).data(c).yFineUB(q,ix) = ar.model(m).data(c).yFineSimu(q,ix);
            q = ar.model(m).data(c).yFineSimu(:,ix) < ar.model(m).data(c).yFineLB(:,ix);
            ar.model(m).data(c).yFineLB(q,ix) = ar.model(m).data(c).yFineSimu(q,ix);
        else
            if(ar.model(m).data(c).has_tExp)
                q = (ar.model(m).data(c).yExpSimu(:,ix) + ar.model(m).data(c).ystdExpSimu(:,ix)) > ar.model(m).data(c).yExpUB(:,ix);
                ar.model(m).data(c).yExpUB(q,ix) = ar.model(m).data(c).yExpSimu(q,ix) + ar.model(m).data(c).ystdExpSimu(q,ix);
                q = (ar.model(m).data(c).yExpSimu(:,ix) - ar.model(m).data(c).ystdExpSimu(:,ix)) < ar.model(m).data(c).yExpLB(:,ix);
                ar.model(m).data(c).yExpLB(q,ix) = ar.model(m).data(c).yExpSimu(q,ix) - ar.model(m).data(c).ystdExpSimu(q,ix);
            end
            q = (ar.model(m).data(c).yFineSimu(:,ix) + ar.model(m).data(c).ystdFineSimu(:,ix)) > ar.model(m).data(c).yFineUB(:,ix);
            ar.model(m).data(c).yFineUB(q,ix) = ar.model(m).data(c).yFineSimu(q,ix) + ar.model(m).data(c).ystdFineSimu(q,ix);
            q = (ar.model(m).data(c).yFineSimu(:,ix) - ar.model(m).data(c).ystdFineSimu(:,ix)) < ar.model(m).data(c).yFineLB(:,ix);
            ar.model(m).data(c).yFineLB(q,ix) = ar.model(m).data(c).yFineSimu(q,ix) - ar.model(m).data(c).ystdFineSimu(q,ix);
        end
            
    catch exception
        fprintf('ERROR for parameter set #%i: %s\n', i, exception.message);
    end
end
 
 
close(hbar)
ar.p = pReset;
ar.qFit = qFit;
arSimu(false, true);
arSimu(false, false);