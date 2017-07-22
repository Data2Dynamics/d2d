% Calculate sensitivities for fitting by finite difference 
%
% arFiniteDifferences(dp)    
%   dp:     parameter variation             [1e-3]


function arFiniteDifferences(dp)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('dp', 'var'))
    dp = 1e-3;
end

if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

pRef = ar.p;
arCalcMerit(false,[]);

ar.sresFD = ar.res' * ones(1,length(pRef)) + 0;
if ( ~isempty( ar.constr ) )
    ar.sconstrFD = ar.constr' * ones(1,length(pRef)) + 0;
end

arWaitbar(0);
for jp=1:length(ar.pLabel)
    arWaitbar(jp, length(ar.pLabel));
    for jm=1:length(ar.model)
        for jc=1:length(ar.model(jm).condition)
            if(jp==1)
                ar.model(jm).condition(jc).suExpSimuFD = zeros(size(ar.model(jm).condition(jc).suExpSimu));
                ar.model(jm).condition(jc).sxExpSimuFD = zeros(size(ar.model(jm).condition(jc).sxExpSimu));
                ar.model(jm).condition(jc).szExpSimuFD = zeros(size(ar.model(jm).condition(jc).szExpSimu));
            end
            qp = ismember(ar.model(jm).condition(jc).p, ar.pLabel{jp}); %R2013a compatible
            if(sum(qp)>0)
                ar.model(jm).condition(jc).suExpSimuFD(:,:,qp) = ar.model(jm).condition(jc).uExpSimu;
                ar.model(jm).condition(jc).sxExpSimuFD(:,:,qp) = ar.model(jm).condition(jc).xExpSimu;
                ar.model(jm).condition(jc).szExpSimuFD(:,:,qp) = ar.model(jm).condition(jc).zExpSimu;
            end
        end
        for jd=1:length(ar.model(jm).data)
            if(jp==1)
                ar.model(jm).data(jd).syExpSimuFD = zeros(size(ar.model(jm).data(jd).syExpSimu));
                ar.model(jm).data(jd).systdExpSimuFD = zeros(size(ar.model(jm).data(jd).systdExpSimu));
                ar.model(jm).data(jd).sresFD = zeros(size(ar.model(jm).data(jd).sres));
                 if( (ar.config.useFitErrorMatrix==0 && ar.config.fiterrors==1) || ...
                        (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(jm,jd)==1) )
                    ar.model(jm).data(jd).sreserrFD = zeros(size(ar.model(jm).data(jd).sres));
                end
            end
            qp = ismember(ar.model(jm).data(jd).p, ar.pLabel{jp}); %R2013a compatible
            if(sum(qp)>0)
                ar.model(jm).data(jd).syExpSimuFD(:,:,qp) = ar.model(jm).data(jd).yExpSimu;
                ar.model(jm).data(jd).systdExpSimuFD(:,:,qp) = ar.model(jm).data(jd).ystdExpSimu;
                ar.model(jm).data(jd).sresFD(:,:,qp) = ar.model(jm).data(jd).res;
                if( (ar.config.useFitErrorMatrix==0 && ar.config.fiterrors==1) || ...
                        (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(jm,jd)==1) )
                    ar.model(jm).data(jd).sreserrFD(:,:,qp) = ar.model(jm).data(jd).reserr;
                end
            end
        end
    end
    
    % perturb p(jp)
    ar.p(jp) = ar.p(jp) + dp;
    try
        arCalcMerit(false,[]);
    catch
        error( 'Failed simulation at %d', jp );
    end
    
    for jm=1:length(ar.model)
        for jc=1:length(ar.model(jm).condition)
            qp = ismember(ar.model(jm).condition(jc).p, ar.pLabel{jp}); %R2013a compatible
            if(sum(qp)>0)
                ar.model(jm).condition(jc).suExpSimuFD(:,:,qp) = ...
                    (ar.model(jm).condition(jc).uExpSimu - ar.model(jm).condition(jc).suExpSimuFD(:,:,qp)) / dp;
                ar.model(jm).condition(jc).sxExpSimuFD(:,:,qp) = ...
                    (ar.model(jm).condition(jc).xExpSimu - ar.model(jm).condition(jc).sxExpSimuFD(:,:,qp)) / dp;
                ar.model(jm).condition(jc).szExpSimuFD(:,:,qp) = ...
                    (ar.model(jm).condition(jc).zExpSimu - ar.model(jm).condition(jc).szExpSimuFD(:,:,qp)) / dp;
            end
        end
        for jd=1:length(ar.model(jm).data)
            qp = ismember(ar.model(jm).data(jd).p, ar.pLabel{jp}); %R2013a compatible
            if(sum(qp)>0)
                ar.model(jm).data(jd).syExpSimuFD(:,:,qp) = ...
                    (ar.model(jm).data(jd).yExpSimu - ar.model(jm).data(jd).syExpSimuFD(:,:,qp)) / dp;
                ar.model(jm).data(jd).systdExpSimuFD(:,:,qp) = ...
                    (ar.model(jm).data(jd).ystdExpSimu - ar.model(jm).data(jd).systdExpSimuFD(:,:,qp)) / dp;
                ar.model(jm).data(jd).sresFD(:,:,qp) = ...
                    (ar.model(jm).data(jd).res - ar.model(jm).data(jd).sresFD(:,:,qp)) / dp;
                if( (ar.config.useFitErrorMatrix==0 && ar.config.fiterrors==1) || ...
                        (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(jm,jd)==1) )
                    ar.model(jm).data(jd).sreserrFD(:,:,qp) = ...
                        (ar.model(jm).data(jd).reserr - ar.model(jm).data(jd).sreserrFD(:,:,qp)) / dp;
                end
            end
        end
    end
    
    ar.sresFD(:,jp) = (ar.res' - ar.sresFD(:,jp)) / dp + 0;
    if ( ~isempty( ar.constr ) )
        ar.sconstrFD(:,jp) = (ar.constr' - ar.sconstrFD(:,jp)) / dp + 0;
    end
    
    % reset p
    ar.p = pRef;
    arCalcMerit(false,[]);
end
arWaitbar(-1);