% arQplot(type, [m], [subset])
% arQplot % alternative usage - see example 7)
% 
% this function sets the model states and observables that are to be plotted.
%
% type      0 all qPlots variables are set to zero, i.e. no plots are opened
%           'y' only plots y, i.e. plots created by arPlotY 
%           'x' only plots x, i.e. plots created by arPlotX
%           'xy' plots x AND y
%           'doseresponse'/'dose-response'/'dr'/'dose' only dose response data
%           'qfit1'/'qfit0' only data that is used / not used for fitting
%           'data_def_name' x and y of this specifc data set matching by ar.model.data.name
% 
% m         which models are plotted. Default behaviour: [1:length(ar.model)]
% 
% subset    subset of plots. The meaning depends on the type of plots:
%           'c': subset == which c [ar.model.condition(c)]
%           'highlighted': subset == d [ar.model.data(d)]
%           in all other cases:
%           plots containing data of data set [ar.model.data(subset)]
% 
% Examples:
% 1) Reset qPlotXs, qPlotYs etc to zero:
% arQplot(0)
% arPlot  % no plot produced
% 
% 2) Plot all data of model 1
% arQplot('y',1)
% arPlot
% 
% 3) Only dose-reponse data:
% arQplot('dose-response')
% 
% 4) Only data, not used for fitting:
% arQplot('qfit0')
% 
% 5) x and y of a specifc data set, indicated by ar.model.data.name
% arQplot('data_def_name')
%
% 6) display information about what is currently plotted in cammand window
% arQplot
%
% see also arPlot

function arQplot(type, m, subset)
global ar

if(~exist('type','var') || isempty(type))
    type = '';
end
if(~exist('m','var') || isempty(m))
    m = 1:length(ar.model);
end

jp = cell(size(ar.model));
if(~exist('subset','var') || isempty(subset))
    for jm=m
        jp{jm} = 1:length(ar.model(jm).qPlotXs);
    end
else
    if(isnumeric(subset))
        subset = Array2Cell(subset*ones(size(ar.model)));
    end
    
    
    for jm=m
        jp{jm} = [];
        for j=1:length(ar.model(jm).plot)
            if ~isempty(intersect(ar.model(jm).plot(j).dLink,subset{jm}))
                jp{jm} = [jp{jm},j];
            end
        end
        jp{jm} = unique(jp{jm});
    end
end


for jm=m
    if isfield(ar.model(jm),'data')
        if isfield(ar.model(jm).data,'name')
            dnames = {ar.model(jm).data.name};
        else
            dnames = cell(0);
        end
        switch(lower(type))
            case 'y'
                ar.model(jm).qPlotXs(:) = 0;
                ar.model(jm).qPlotYs(jp{jm}) = 1;
                ar.model(jm).qPlotVs(:) = 0;
            case 'v'
                ar.model(jm).qPlotXs(:) = 0;
                ar.model(jm).qPlotYs(:) = 0;
                ar.model(jm).qPlotVs(jp{jm}) = 1;
            case 'x'
                ar.model(jm).qPlotXs(jp{jm}) = 1;
                ar.model(jm).qPlotYs(:) = 0;
                ar.model(jm).qPlotVs(:) = 0;
            case {'xy','yx'}
                ar.model(jm).qPlotXs(jp{jm}) = 1;
                ar.model(jm).qPlotYs(jp{jm}) = 1;
                ar.model(jm).qPlotVs(:) = 0;
                
            case {'xv','vx'}
                ar.model(jm).qPlotXs(jp{jm}) = 1;
                ar.model(jm).qPlotYs(jp{jm}) = 0;
                ar.model(jm).qPlotVs(:) = 1;
                
            case {'yv','vy'}
                ar.model(jm).qPlotXs(jp{jm}) = 0;
                ar.model(jm).qPlotYs(jp{jm}) = 1;
                ar.model(jm).qPlotVs(:) = 1;
                
            case {'vxy','all','xyv','yxv','yvx','xvy','vyx','on',1,true}
                ar.model(jm).qPlotCs(:) = 0;
                ar.model(jm).qPlotXs(jp{jm}) = 1;
                ar.model(jm).qPlotYs(jp{jm}) = 1;
                ar.model(jm).qPlotVs(jp{jm}) = 1;
                
            case {0,false,'off'}
                ar.model(jm).qPlotXs(jp{jm}) = 0;
                ar.model(jm).qPlotYs(jp{jm}) = 0;
                ar.model(jm).qPlotVs(jp{jm}) = 0;
                
            case {'doseresponse','dr','dose','dose-response'}
                ar.model(jm).qPlotYs(:) = 0;
                for p=jp{jm}
                    if(ar.model(jm).plot(p).doseresponse ==1)
                        ar.model(jm).qPlotYs(p) = 1;
                    end
                end
                
            case {'hl','highlighted','highlight'}
                ar.model(jm).qPlotYs(:) = 0;
                for p=jp{jm}
                    for subset=ar.model(jm).plot(p).dLink
                        if(sum(ar.model(jm).data(subset).highlight(:))>0)
                            ar.model(jm).qPlotYs(p) = 1;
                        end
                    end
                end
                
            case  'qfit0'
                ar.model(jm).qPlotYs(:) = 0;
                for p=jp{jm}
                    for subset=ar.model(jm).plot(p).dLink
                        if(sum(ar.model(jm).data(subset).qFit==0)>0)
                            ar.model(jm).qPlotYs(p) = 1;
                        end
                    end
                end
                
            case  'qfit1'
                ar.model(jm).qPlotYs(:) = 0;
                for p=jp{jm}
                    for subset=ar.model(jm).plot(p).dLink
                        if(sum(ar.model(jm).data(subset).qFit==1)>0)
                            ar.model(jm).qPlotYs(p) = 1;
                        end
                    end
                end
                
            case dnames
                indd = strmatch(type,dnames,'exact');
                
                for p=jp{jm}
                    ind = intersect(indd,ar.model(jm).plot(p).dLink);
                    if(~isempty(ind))
                        ar.model(jm).qPlotXs(p) = 1;
                        ar.model(jm).qPlotYs(p) = 1;
                    end
                end
                
            otherwise
                error('')
        end
    end
end

if(isempty(type))
    for jm=1:length(ar.model)
        disp(['ar.model(',num2str(jm),'):'])
        disp(['find(qPlotVs) = ',sprintf('%i, ',find(ar.model(jm).qPlotVs))]);
        disp(['find(qPlotXs) = ',sprintf('%i, ',find(ar.model(jm).qPlotXs))]);
        disp(['find(qPlotYs) = ',sprintf('%i, ',find(ar.model(jm).qPlotYs))]);
    end
end


function c=Array2Cell(a)
c = arrayfun(@(x){x},a);
