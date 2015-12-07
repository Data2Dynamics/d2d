function varargout = arPrintChi2(doYs, onesided)

if(~exist('doYs','var'))
    doYs = false;
end
if(~exist('onesided','var'))
    onesided = false;
end

if(doYs)
    [a,b,c] = tmpforYs(nargout==0, onesided);
    if(nargout>0)
        varargout{1} = a;
        varargout{2} = b;
        varargout{3} = c;
    end
    return
end

global ar

if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

a = [];
b = {};
for jm=1:length(ar.model)
    for jplot = 1:length(ar.model(jm).plot)
        chi2 = 0;
        chi2err = 0;
        ndata = 0;
        ndata_err = 0;
        for jd = ar.model(jm).plot(jplot).dLink
            chi2 = chi2 + sum(ar.model(jm).data(jd).chi2(ar.model(jm).data(jd).qFit==1));
            ndata = ndata + sum(ar.model(jm).data(jd).ndata(ar.model(jm).data(jd).qFit==1));
            if(ar.config.useFitErrorMatrix == 0 && ar.config.fiterrors == 1)
                chi2err = chi2err + sum(ar.model(jm).data(jd).chi2err(ar.model(jm).data(jd).qFit==1));
            elseif(ar.config.useFitErrorMatrix == 1 && ar.config.fiterrors_matrix(jm,jd)==1)
                chi2err = chi2err + sum(ar.model(jm).data(jd).chi2err(ar.model(jm).data(jd).qFit==1));
                ndata_err = ndata_err + sum(ar.model(jm).data(jd).ndata(ar.model(jm).data(jd).qFit==1));
            end
        end
        chi2fit = chi2 + chi2err;
        a(end+1) = chi2fit;
        b{end+1} = [ar.model(jm).name '_' ar.model(jm).plot(jplot).name];
        
        if(nargout==0)
            if(ar.config.useFitErrorMatrix == 0 && ar.config.fiterrors == 1)
                fprintf('m=%20s, d=%50s :\t -2*log(L) = %8.2g, %4i data points\n', ...
                    ar.model(jm).name, ar.model(jm).plot(jplot).name, ...
                    2*ndata*log(sqrt(2*pi)) + chi2fit, ndata);
            elseif(ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(jm,ar.model(jm).plot(jplot).dLink(1))==1)
                fprintf('m=%20s, d=%50s :\t -2*log(L) = %8.2g, %4i data points\n', ...
                    ar.model(jm).name, ar.model(jm).plot(jplot).name, ...
                    2*ndata_err*log(sqrt(2*pi)) + chi2fit, ndata);
            else
                fprintf('m=%20s, d=%50s :\t chi^2 = %8.2g, %4i data points\n', ...
                    ar.model(jm).name, ar.model(jm).plot(jplot).name, ...
                    chi2, ndata);
            end
        end
    end
end

if(nargout>0)
    varargout{1} = a;
    varargout{2} = b;
end

function [a,b,c] = tmpforYs(dodisp, onesided)

global ar

a = [];
b = {};
c = [];
for jm=1:length(ar.model)
    for jplot = 1:length(ar.model(jm).plot)
        for jy = 1:ar.model(jm).plot(jplot).ny
            chi2 = 0;
            chi2err = 0;
            ndata = 0;
            ndata_err = 0;
            res = [];
            for jd = ar.model(jm).plot(jplot).dLink
                if(ar.model(jm).data(jd).qFit(jy)==1)
                    chi2 = chi2 + sum(ar.model(jm).data(jd).chi2(jy));
                    ndata = ndata + sum(ar.model(jm).data(jd).ndata(jy));
                    res = [res; ar.model(jm).data(jd).res(:,jy)];
                    if(ar.config.useFitErrorMatrix == 0 && ar.config.fiterrors == 1)
                        chi2err = chi2err + sum(ar.model(jm).data(jd).chi2err(ar.model(jm).data(jd).qFit==1));
                    elseif(ar.config.useFitErrorMatrix == 1 && ar.config.fiterrors_matrix(jm,jd)==1)
                        chi2err = chi2err + sum(ar.model(jm).data(jd).chi2err(ar.model(jm).data(jd).qFit==1));
                        ndata_err = ndata_err + sum(ar.model(jm).data(jd).ndata(ar.model(jm).data(jd).qFit==1));
                    end
                end
            end
            chi2fit = chi2 + chi2err;
            a(end+1) = chi2fit;
            b{end+1} = [ar.model(jm).name '_' ar.model(jm).plot(jplot).name '_' ar.model(jm).data(jd).y{jy}];
            res_var = sum(res.^2) / ndata;
            c(end+1) = res_var;
            if(~onesided)
                q = chi2inv([0.05 0.95],ndata)/ndata;
                if(q(1)<res_var && res_var<q(2))
                    flag = '';
                else
                    flag = '(***)';
                end
            else
                q = chi2inv(0.95,ndata)/ndata;
                if(res_var < q)
                    flag = '';
                else
                    flag = '(***)';
                end
            end

            
            if(dodisp)
                 if(ar.config.useFitErrorMatrix == 0 && ar.config.fiterrors == 1)
                    fprintf('m=%20s, d=%50s, y=%20s :\t -2*log(L) = %8.2g, %4i data points, var(res) = %f %s\n', ...
                        ar.model(jm).name, ar.model(jm).plot(jplot).name, ar.model(jm).data(jd).y{jy}, ...
                        2*ndata*log(sqrt(2*pi)) + chi2fit, ndata, res_var, flag);
                elseif(ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(jm,ar.model(jm).plot(jplot).dLink(1))==1)
                    fprintf('m=%20s, d=%50s, y=%20s :\t -2*log(L) = %8.2g, %4i data points, var(res) = %f %s\n', ...
                        ar.model(jm).name, ar.model(jm).plot(jplot).name, ar.model(jm).data(jd).y{jy}, ...
                        2*ndata_err*log(sqrt(2*pi)) + chi2fit, ndata, res_var, flag);
                else
                    fprintf('m=%20s, d=%50s, y=%20s :\t chi^2 = %8.2g, %4i data points, var(res) = %f %s\n', ...
                        ar.model(jm).name, ar.model(jm).plot(jplot).name, ar.model(jm).data(jd).y{jy}, ...
                        chi2, ndata, res_var, flag);
                end
            end
        end
    end
end

if(dodisp)
    [~, ri] = sort(c);
    for j=1:length(c);
        fprintf('var(res) = %f\t %60s\n', c(ri(j)), b{ri(j)});
    end
end


