% arCalcRes(sensi)
% 
%   This function performs calculation of residuals (res and reserr) for
%   the individual data sets including their derivatives sres, sreserr
% 
%   In this function, there is the major evaluation of ar.config.fiterrors
%   i.e. here the decision is made whether reserr has to be calculated and
%   whether experimental errors or errors from the error model are taken.
% 
%   The function implements this step as previously done in arSimuCalc.c
%   The respective C-code is still in this function as comment.
% 
%   This function also uses the additive constant ar.config.add_c=50 which is
%   required for using lsqnonlin in case of fitting errors.
% 
%       sensi   should sensitivities sres, sreserr be calculated?
%               Default: sensi=1


function arCalcRes(sensi)
if ~exist('sensi','var') || isempty(sensi)
    sensi = 1;
end

global ar

if(~isfield(ar.config,'useFitErrorCorrection'))
    ar.config.useFitErrorCorrection = true;
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end
if ar.config.useFitErrorMatrix==1
    ar.ndata_err = 0;
end

ar.ndata_res = 0;

% correction for error fitting
for jm = 1:length(ar.model)
    if(isfield(ar.model(jm), 'data'))
        nd = length(ar.model(jm).data);
        for jd = 1:nd
            if(ar.model(jm).data(jd).has_yExp)
                ar.ndata_res = ar.ndata_res + sum(ar.model(jm).data(jd).ndata(ar.model(jm).data(jd).qFit==1));
                if(ar.config.useFitErrorMatrix == 1 && ar.config.fiterrors_matrix(jm,jd) == 1)
                    ar.ndata_err = ar.ndata_err + sum(ar.model(jm).data(jd).ndata(ar.model(jm).data(jd).qFit==1));
                end
            end
        end
    end
end

% if error parameter fitted or fixed, then useFitErrorCorrection is
% evaluated:
if  ar.ndata_res>0 && (ar.config.fiterrors==1 || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)>0) && ar.config.useFitErrorCorrection  )
    if(ar.ndata_res -sum(ar.qError~=1 & ar.qFit==1) < sum(ar.qError~=1 & ar.qFit==1))
        ar.config.fiterrors_correction = 1;
        if(~ar.config.fiterrors_correction_warning)
            warning('ar.config.fiterrors_correction_warning : turning off bias correction, not enough data'); %#ok<WNTAG>
            ar.config.fiterrors_correction_warning = true;
        end
    else
        ar.config.fiterrors_correction = ar.ndata_res/(ar.ndata_res-sum(ar.qError~=1 & ar.qFit==1));
        ar.config.fiterrors_correction_warning = false;
    end
else
    ar.config.fiterrors_correction = 1;
end

%% add user-defined residual(s)
% If you want to modify the objective function by adding residuals, use a
% user-defined function specified as ar.config.user_residual_fun
% If this field is available, then the function is called here:
ar.res_user = struct;
if isfield(ar.config,'user_residual_fun') && ~isempty(ar.config.user_residual_fun)
    [ar.res_user.res,ar.res_user.sres,ar.res_user.type] = feval(ar.config.user_residual_fun); % this function can implement additional residuals added to ar.res
    if length(ar.res_user.res)~=size(ar.res_user.sres,1)
        error('length(ar.res_user.res)~=size(ar.res_user.sres,1)')
    end
else
    ar.res_user.res = [];
    ar.res_user.sres = [];
    ar.res_user.type = [];
end


fiterrors_correction_factor = ar.config.fiterrors_correction;

for m=1:length(ar.model)
    if ( isfield( ar.model(m), 'data' ) )
        for d=1:length(ar.model(m).data)

           %% this is THE point in D2D where ar.config.fiterrors enters:
            % ar.config.fiterrors == 0: use exp. errors where available and error
            % model otherwise.
            % Fitting of error parameters is exclusively controlled by ar.qFit 
            if ar.config.fiterrors==0
                ystd = ar.model(m).data(d).yExpStd;

                noSD = isnan(ystd);
                ystd(noSD) = ar.model(m).data(d).ystdExpSimu(noSD); % not available => use error model

                systd = zeros(size(ar.model(m).data(d).systdExpSimu)); % if experimental error available => no parameter dependency            
                for ip=1:size(systd,3)
                    tmp = zeros(size(ystd)); 
                    syExpSimu = ar.model(m).data(d).systdExpSimu(:,:,ip);
                    tmp(noSD) = syExpSimu(noSD);  
                    systd(:,:,ip) = tmp;               
                end

            % ar.config.fiterrors==-1: only use exp. errors
            elseif ar.config.fiterrors == -1  % ensure that only experimental errors are used
                ystd = ar.model(m).data(d).yExpStd;
                systd = zeros(size(ar.model(m).data(d).systdExpSimu));


            % ar.config.fiterrors == 1: only use exp. error model (and omit exp.
            % errors)
            elseif ar.config.fiterrors == 1 % only error model is used
                ystd = ar.model(m).data(d).ystdExpSimu;
                systd = ar.model(m).data(d).systdExpSimu;               

            else
                error('ar.config.fiterrors = %f is not yet implemented',ar.config.fiterrors);
            end


            [ar.model(m).data(d).res, ar.model(m).data(d).chi2]         = fres(ar.model(m).data(d).yExp, ar.model(m).data(d).yExpSimu, ystd, fiterrors_correction_factor);
            if isempty(ar.model(m).data(d).res)
                ar.model(m).data(d).chi2 = zeros(1,length(ar.model(m).data(d).fy));
            end

            if (ar.config.fiterrors == 1) || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)>0 ) % error residuals are only !=0 if errors are fitted:
                [ar.model(m).data(d).reserr, ar.model(m).data(d).chi2err]   = fres_error(ystd, ar.config.add_c);
                if isempty(ar.model(m).data(d).reserr)
                    ar.model(m).data(d).chi2err = zeros(1,length(ar.model(m).data(d).fy));
                end
                    
            else
                ar.model(m).data(d).reserr = zeros(size(ar.model(m).data(d).res));
                ar.model(m).data(d).chi2err = zeros(size(ar.model(m).data(d).chi2));
            end

            if (sensi == 1) && ar.model(m).data(d).has_yExp 
                ar.model(m).data(d).sres = fsres(ar.model(m).data(d).syExpSimu, ar.model(m).data(d).yExp, ystd, fiterrors_correction_factor);
                if (ar.config.fiterrors == 1) || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)>0 )
                    [ar.model(m).data(d).sres, ar.model(m).data(d).sreserr] = fsres_error(ar.model(m).data(d).yExp, ar.model(m).data(d).res, ar.model(m).data(d).reserr, ar.model(m).data(d).sres, ystd, systd);
                end
            else
                ar.model(m).data(d).sres(:) = NaN;
                ar.model(m).data(d).sreserr(:) = NaN;
            end



            %         /* log trafo of parameters */

            if sensi && ar.model(m).data(d).has_yExp
                for ip=1:size(ar.model(m).data(d).sres,3)
                    if ar.model(m).data(d).qLog10(ip) > 0.5
                        ar.model(m).data(d).sres(:,:,ip)    = ar.model(m).data(d).sres(:,:,ip).*ar.model(m).data(d).pNum(ip) * log(10);
                        if (ar.config.fiterrors == 1)  || (ar.config.fiterrors==0 && sum(ar.qFit(ar.qError==1)<2)>0 )
                            ar.model(m).data(d).sreserr(:,:,ip) = ar.model(m).data(d).sreserr(:,:,ip).*ar.model(m).data(d).pNum(ip) * log(10);
                        end
                        %             for (ip=0; ip < np; ip++) {
                        %                 if (qlogp[ip] > 0.5) {
                        %                      for (iy=0; iy<ny; iy++) {
                        %                          sres[it + (iy*nt) + (ip*nt*ny)] *= p[ip] * log(10.0);
                        %                          if( (useFitErrorMatrix == 0 && fiterrors == 1) || (useFitErrorMatrix == 1 && fiterrors_matrix[id*nrows_fiterrors_matrix+im] == 1) ) {
                        %                              sreserr[it + (iy*nt) + (ip*nt*ny)] *= p[ip] * log(10.0);
                        %                          }
                        %                      }
                    end
                end
            end
        end
    end
end

end

%  /* standard least squares */
function [res,chi2] = fres(yexp, y, ystd, fiterrors_correction_factor)
    res = (yexp-y)./ystd * sqrt(fiterrors_correction_factor);
    res(isnan(res)) = 0;
    chi2 = sum(res.^2,1);
end
%  void fres(int nt, int ny, int it, double *res, double *y, double *yexp, double *ystd, double *chi2, double fiterrors_correction_factor) {
%      int iy;
%      
%      for(iy=0; iy<ny; iy++){
%          res[it + (iy*nt)] = (yexp[it + (iy*nt)] - y[it + (iy*nt)]) / ystd[it + (iy*nt)] * sqrt(fiterrors_correction_factor);
%          /* in case of missing data (nan) */
%          if(mxIsNaN(yexp[it + (iy*nt)])) {
%              res[it + (iy*nt)] = 0.0;
%              y[it + (iy*nt)] = yexp[it + (iy*nt)];
%              ystd[it + (iy*nt)] = yexp[it + (iy*nt)];
%          }
%          /* in case of Inf data after log10(0) */
%          if(mxIsInf(yexp[it + (iy*nt)])) {
%              res[it + (iy*nt)] = 0.0;
%          }
%          chi2[iy] += pow(res[it + (iy*nt)], 2);
%      }
%  }

function sres = fsres(sy, yexp, ystd, fiterrors_correction_factor)
sres = NaN(size(sy));
    for ip=1:size(sy,3)
        sres(:,:,ip) = - sy(:,:,ip) ./ ystd * sqrt(fiterrors_correction_factor);
        tmp = sres(:,:,ip);
        tmp(isnan(yexp)) = 0;
        tmp(isinf(yexp)) = 0;
        sres(:,:,ip) = tmp;
    end
end
%  void fsres(int nt, int ny, int np, int it, double *sres, double *sy, double *yexp, double *ystd, double fiterrors_correction_factor) {
%      int iy, ip;
%      
%      for(iy=0; iy<ny; iy++){
%          for(ip=0; ip<np; ip++){
%              sres[it + (iy*nt) + (ip*nt*ny)] = - sy[it + (iy*nt) + (ip*nt*ny)] / ystd[it + (iy*nt)] * sqrt(fiterrors_correction_factor);
%              /* in case of missing data (nan) */
%              if(mxIsNaN(yexp[it + (iy*nt)])) {
%                  sres[it + (iy*nt) + (ip*nt*ny)] = 0.0;
%              }
%              /* in case of Inf data after log10(0) */
%              if(mxIsInf(yexp[it + (iy*nt)])) {
%                  sres[it + (iy*nt) + (ip*nt*ny)] = 0.0;
%              }
%          }
%      }
%  }
%  
%  /* least squares for error model fitting */
function [reserr,chi2err] = fres_error(ystd, add_c)

    reserr = 2.0*log(ystd) + add_c;    
    reserr(isnan(ystd)) = 0;
    
    if(sum(reserr(:) < 0)>0)
        error('arCalcRes/fres_error: error residual too small.');
    else 
        reserr = sqrt(reserr);
        chi2err = sum(reserr.^2,1) - add_c*sum(abs(reserr)>0,1);  
    end    

%     reserr
%     chi2err
%     
%      add_c = 50.0;
%      reserr = zeros(size(ystd));
%      chi2err = zeros(size(ystd(1,:)));
%      
%      for it=1:size(yexp,1)
%          for iy=1:size(ystd,2)
%              reserr(it,iy) = 2.0*log(ystd(it,iy));
%              if(isnan(yexp(it,iy)))
%                  reserr(it,iy) = 0.0;
%                  y(it,iy) = yexp(it,iy);
%                  ystd(it,iy) = yexp(it,iy);
%              else
%                  reserr(it,iy) = reserr(it,iy) + add_c;
%                  
%                  if(reserr(it,iy) < 0)
%                      fprintf('ERROR error model < 1e-10 not allowed\n');
%                      return;
%                  end
%                  reserr(it,iy) = sqrt(reserr(it,iy));
%                  chi2err(iy) = chi2err(iy) + reserr(it,iy).^2 - add_c;
%              end
%          end
%      end
%      
%      reserr
%      chi2err
end
    %  void fres_error(int nt, int ny, int it, double *reserr, double *res, double *y, double *yexp, double *ystd, double *chi2err) {
%      int iy;
%      
%      double add_c = 50.0;
%      
%      for(iy=0; iy<ny; iy++){
%          reserr[it + (iy*nt)] = 2.0*log(ystd[it + (iy*nt)]);
%          if(mxIsNaN(yexp[it + (iy*nt)])) {
%              reserr[it + (iy*nt)] = 0.0;
%              y[it + (iy*nt)] = yexp[it + (iy*nt)];
%              ystd[it + (iy*nt)] = yexp[it + (iy*nt)];
%          } else {
%              reserr[it + (iy*nt)] += add_c;
%              /* 2*log(ystd) + add_c > 0 */
%              if(reserr[it + (iy*nt)] < 0) {
%                  printf("ERROR error model < 1e-10 not allowed\n");
%                  return;
%              }
%              reserr[it + (iy*nt)] = sqrt(reserr[it + (iy*nt)]);
%              chi2err[iy] += pow(reserr[it + (iy*nt)], 2) - add_c;
%          }
%      }
%  }
function [sres,sreserr] = fsres_error(yexp, res, reserr, sres, ystd, systd)
sreserr = NaN(size(sres));
for ip=1:size(sres,3)
    sres(:,:,ip) = sres(:,:,ip) - (systd(:,:,ip) .* res ./ ystd);
    sreserr(:,:,ip) = systd(:,:,ip) ./ (reserr.*ystd);
    
    tmp = sres(:,:,ip);
    tmp(isnan(yexp)) = 0;
    sres(:,:,ip) = tmp;
    
    tmp = sreserr(:,:,ip);
    tmp(isnan(yexp)) = 0;
    sreserr(:,:,ip) = tmp;
end
end

%  void fsres_error(int nt, int ny, int np, int it, double *sres, double *sreserr, double *sy, double *systd, double *y, double *yexp, double *ystd, double *res, double *reserr) {
%      int iy, ip;
%      
%      for(iy=0; iy<ny; iy++){
%          for(ip=0; ip<np; ip++){
%              sres[it + (iy*nt) + (ip*nt*ny)] -= systd[it + (iy*nt) + (ip*nt*ny)] * res[it + (iy*nt)] / ystd[it + (iy*nt)];
%              sreserr[it + (iy*nt) + (ip*nt*ny)] = systd[it + (iy*nt) + (ip*nt*ny)] / (reserr[it + (iy*nt)] * ystd[it + (iy*nt)]);
%              if(mxIsNaN(yexp[it + (iy*nt)])) {
%                  sres[it + (iy*nt) + (ip*nt*ny)] = 0.0;
%                  sreserr[it + (iy*nt) + (ip*nt*ny)] = 0.0;
%              }
%          }
%      }
%  }
%  

