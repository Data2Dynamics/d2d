% plot confidence bands of trajectories of ple
%
% arPLEConfidence(jks, n, do_vali)
%
% jks               par response                                    [all]
% n                 trajectories per parameter                      [10] (0=all)
% do_vali           compute validation confidence bands             [false]
% do_refit_obs      refit observation parameters                    [false]
function arPLEConfidence(jks, n, do_vali, do_refit_obs)

global ar
global pleGlobals

if(isempty(ar))
	error('please initialize by arInit')
end
if(isempty(pleGlobals))
	error('PLE ERROR: please initialize')
end

if(~exist('jks','var') || isempty(jks))
	jks = find(ar.qFit==1);
end
if(~exist('n','var') || isempty(n))
	n = 10;
end
if(~exist('do_vali','var'))
	do_vali = false;
end
if(~exist('do_refit_obs','var'))
	do_refit_obs = false;
end

jks = find(ismember(pleGlobals.p_labels, ar.pLabel(jks))); %R2013a compatible

for jm=1:length(ar.model)
	for jc = 1:length(ar.model(jm).condition)
		ar.model(jm).condition(jc).uExpUB = ar.model(jm).condition(jc).uExpSimu;
		ar.model(jm).condition(jc).uExpLB = ar.model(jm).condition(jc).uExpSimu;
		ar.model(jm).condition(jc).uFineUB = ar.model(jm).condition(jc).uFineSimu;
		ar.model(jm).condition(jc).uFineLB = ar.model(jm).condition(jc).uFineSimu;

        ar.model(jm).condition(jc).vExpUB = ar.model(jm).condition(jc).vExpSimu;
		ar.model(jm).condition(jc).vExpLB = ar.model(jm).condition(jc).vExpSimu;
		ar.model(jm).condition(jc).vFineUB = ar.model(jm).condition(jc).vFineSimu;
		ar.model(jm).condition(jc).vFineLB = ar.model(jm).condition(jc).vFineSimu;
        
		ar.model(jm).condition(jc).xExpUB = ar.model(jm).condition(jc).xExpSimu;
		ar.model(jm).condition(jc).xExpLB = ar.model(jm).condition(jc).xExpSimu;
		ar.model(jm).condition(jc).xFineUB = ar.model(jm).condition(jc).xFineSimu;
		ar.model(jm).condition(jc).xFineLB = ar.model(jm).condition(jc).xFineSimu;
	end
	if(isfield(ar.model(jm), 'data'))
		for jd = 1:length(ar.model(jm).data)
            if(ar.model(jm).data(jd).has_tExp)
                ar.model(jm).data(jd).yExpUB = ar.model(jm).data(jd).yExpSimu;
                ar.model(jm).data(jd).yExpLB = ar.model(jm).data(jd).yExpSimu;
            end
			ar.model(jm).data(jd).yFineUB = ar.model(jm).data(jd).yFineSimu;
			ar.model(jm).data(jd).yFineLB = ar.model(jm).data(jd).yFineSimu;
		end
	end
end

% for reseting
pReset = ar.p;

% collect
hbar = waitbar(0, 'Please wait...');

for j = jks
	if(j<=length(pleGlobals.chi2s) && ~isempty(pleGlobals.ps{j}))
		% all non nans
		chi2stmp = pleGlobals.chi2s{j}(~isnan(pleGlobals.chi2s{j}));
		pstmp = pleGlobals.ps{j}(~isnan(pleGlobals.chi2s{j}),:);

%         % cut out bounds
%         jk = j;
% 		qCloseToUB = pstmp > ones(length(chi2stmp),1) * (pleGlobals.ub - pleGlobals.dist_thres) & ...
% 			ones(length(chi2stmp),1) * pleGlobals.q_fit==1;
% 		qCloseToLB = pstmp < ones(length(chi2stmp),1) * (pleGlobals.lb + pleGlobals.dist_thres) & ...
% 			ones(length(chi2stmp),1) * pleGlobals.q_fit==1;
% 
% 		qhitbound = false(size(pstmp));
% 		qhitbound(:,pleGlobals.q_fit==1) = pleGlobals.gradient{jk}(~isnan(pleGlobals.chi2s{j}),pleGlobals.q_fit==1) > pleGlobals.grad_thres & qCloseToLB(:,pleGlobals.q_fit==1) | ...
% 			pleGlobals.gradient{jk}(~isnan(pleGlobals.chi2s{j}),pleGlobals.q_fit==1) < -pleGlobals.grad_thres & qCloseToUB(:,pleGlobals.q_fit==1);
% 
% 		chi2stmp = chi2stmp(sum(qhitbound,2)==0);
% 		pstmp = pstmp(sum(qhitbound,2)==0,:);

		% cut below point-wise threshold
		qbelow = chi2stmp < pleGlobals.chi2 + pleGlobals.dchi2_point;
		chi2stmp = chi2stmp(qbelow);
		pstmp = pstmp(qbelow,:);

		% make n subset
		if(n>0 && length(chi2stmp) > n)
			iselectn = floor(linspace(1, length(chi2stmp), n));
			pstmp = pstmp(iselectn,:);
		end

		for jp = 1:size(pstmp,1)
			hbar = waitbar(((j-1)*size(pstmp,1) + jp)/(size(pstmp,1)*length(jks)), hbar, 'Please wait...');
            ar.p(ismember(ar.pLabel, pleGlobals.p_labels)) = pstmp(jp,ismember(pleGlobals.p_labels, ar.pLabel)); %R2013a compatible
			try
                if(do_refit_obs)
                    arFitObs(true);
                end
				arSimu(false, true);
				arSimu(false, false);

				for jm=1:length(ar.model)
					for jc = 1:length(ar.model(jm).condition)
						q = ar.model(jm).condition(jc).uExpSimu > ar.model(jm).condition(jc).uExpUB;
						ar.model(jm).condition(jc).uExpUB(q) = ar.model(jm).condition(jc).uExpSimu(q);
						q = ar.model(jm).condition(jc).uExpSimu < ar.model(jm).condition(jc).uExpLB;
						ar.model(jm).condition(jc).uExpLB(q) = ar.model(jm).condition(jc).uExpSimu(q);
						q = ar.model(jm).condition(jc).uFineSimu > ar.model(jm).condition(jc).uFineUB;
						ar.model(jm).condition(jc).uFineUB(q) = ar.model(jm).condition(jc).uFineSimu(q);
						q = ar.model(jm).condition(jc).uFineSimu < ar.model(jm).condition(jc).uFineLB;
						ar.model(jm).condition(jc).uFineLB(q) = ar.model(jm).condition(jc).uFineSimu(q);
                        
						q = ar.model(jm).condition(jc).vExpSimu > ar.model(jm).condition(jc).vExpUB;
						ar.model(jm).condition(jc).vExpUB(q) = ar.model(jm).condition(jc).vExpSimu(q);
						q = ar.model(jm).condition(jc).vExpSimu < ar.model(jm).condition(jc).vExpLB;
						ar.model(jm).condition(jc).vExpLB(q) = ar.model(jm).condition(jc).vExpSimu(q);
						q = ar.model(jm).condition(jc).vFineSimu > ar.model(jm).condition(jc).vFineUB;
						ar.model(jm).condition(jc).vFineUB(q) = ar.model(jm).condition(jc).vFineSimu(q);
						q = ar.model(jm).condition(jc).vFineSimu < ar.model(jm).condition(jc).vFineLB;
						ar.model(jm).condition(jc).vFineLB(q) = ar.model(jm).condition(jc).vFineSimu(q);

						q = ar.model(jm).condition(jc).xExpSimu > ar.model(jm).condition(jc).xExpUB;
						ar.model(jm).condition(jc).xExpUB(q) = ar.model(jm).condition(jc).xExpSimu(q);
						q = ar.model(jm).condition(jc).xExpSimu < ar.model(jm).condition(jc).xExpLB;
						ar.model(jm).condition(jc).xExpLB(q) = ar.model(jm).condition(jc).xExpSimu(q);
						q = ar.model(jm).condition(jc).xFineSimu > ar.model(jm).condition(jc).xFineUB;
						ar.model(jm).condition(jc).xFineUB(q) = ar.model(jm).condition(jc).xFineSimu(q);
						q = ar.model(jm).condition(jc).xFineSimu < ar.model(jm).condition(jc).xFineLB;
						ar.model(jm).condition(jc).xFineLB(q) = ar.model(jm).condition(jc).xFineSimu(q);
					end
					if(isfield(ar.model(jm), 'data'))
						for jd = 1:length(ar.model(jm).data)
                            if(~do_vali)
                                if(ar.model(jm).data(jd).has_tExp)
                                    q = ar.model(jm).data(jd).yExpSimu > ar.model(jm).data(jd).yExpUB;
                                    ar.model(jm).data(jd).yExpUB(q) = ar.model(jm).data(jd).yExpSimu(q);
                                    q = ar.model(jm).data(jd).yExpSimu < ar.model(jm).data(jd).yExpLB;
                                    ar.model(jm).data(jd).yExpLB(q) = ar.model(jm).data(jd).yExpSimu(q);
                                end
                                q = ar.model(jm).data(jd).yFineSimu > ar.model(jm).data(jd).yFineUB;
                                ar.model(jm).data(jd).yFineUB(q) = ar.model(jm).data(jd).yFineSimu(q);
                                q = ar.model(jm).data(jd).yFineSimu < ar.model(jm).data(jd).yFineLB;
                                ar.model(jm).data(jd).yFineLB(q) = ar.model(jm).data(jd).yFineSimu(q);
                            else
                                if(ar.model(jm).data(jd).has_tExp)
                                    q = (ar.model(jm).data(jd).yExpSimu + ar.model(jm).data(jd).ystdExpSimu) > ar.model(jm).data(jd).yExpUB;
                                    ar.model(jm).data(jd).yExpUB(q) = ar.model(jm).data(jd).yExpSimu(q) + ar.model(jm).data(jd).ystdExpSimu(q);
                                    q = (ar.model(jm).data(jd).yExpSimu - ar.model(jm).data(jd).ystdExpSimu) < ar.model(jm).data(jd).yExpLB;
                                    ar.model(jm).data(jd).yExpLB(q) = ar.model(jm).data(jd).yExpSimu(q) - ar.model(jm).data(jd).ystdExpSimu(q);
                                end
                                q = (ar.model(jm).data(jd).yFineSimu + ar.model(jm).data(jd).ystdFineSimu) > ar.model(jm).data(jd).yFineUB;
                                ar.model(jm).data(jd).yFineUB(q) = ar.model(jm).data(jd).yFineSimu(q) + ar.model(jm).data(jd).ystdFineSimu(q);
                                q = (ar.model(jm).data(jd).yFineSimu - ar.model(jm).data(jd).ystdFineSimu) < ar.model(jm).data(jd).yFineLB;
                                ar.model(jm).data(jd).yFineLB(q) = ar.model(jm).data(jd).yFineSimu(q) - ar.model(jm).data(jd).ystdFineSimu(q);
                            end
						end
					end
				end

            catch exception
				fprintf('ERROR for parameter set #%i: %s\n', jp, exception.message);
			end
		end
	end
end

close(hbar)
ar.p = pReset;
arSimu(false, true);
arSimu(false, false);
