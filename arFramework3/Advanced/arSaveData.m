% arSaveData(m, jplot)
%
% Save data sets corresponding to the plot with index jplot to csv file
% into the subfolder /Data
%
%   m      [1:length(ar.model)]           position of model			
%   jplot  [1:length(ar.model(m).plot)]   position of plot
%
% The data file is saved to ['./Data/' ar.model(m).plot(jplot).name '.csv']

function arSaveData(m, jplot)

global ar

if(isempty(ar))
	error('please initialize by arInit')
end

if(~exist('Data', 'dir'))
	mkdir('Data')
end

if(~exist('m','var'))
    for jm=1:length(ar.model)
        for jplot=1:length(ar.model(jm).plot)
            arSaveData(jm, jplot)
        end
    end
    return
end
if(~exist('jplot','var'))
    for jplot = 1:length(ar.model(m).plot)
        arSaveData(m, jplot)
    end
    return
end
if(~isfield(ar.config,'useFitErrorMatrix'))
    ar.config.useFitErrorMatrix = false;
end

fname = ['./Data/' ar.model(m).plot(jplot).name '.csv'];
fprintf('saving data from model #%i, plot #%i, %s...', m, jplot, fname);

fid = fopen(fname, 'w');

jd = ar.model(m).plot(jplot).dLink(1);

% save to file
fprintf(fid, '"%s",', ar.model(m).data(jd).tUnits{3});

% conditions
if(~isempty(ar.model(m).data(jd).condition))
    for jp = 1:length(ar.model(m).data(jd).condition)
        fprintf(fid, '"%s",', ar.model(m).data(jd).condition(jp).parameter);
    end
end

% y headers
fprintf(fid, '"%s",', ar.model(m).data(jd).y{:});

% ystd headers
if( (ar.config.useFitErrorMatrix == 0 && ar.config.fiterrors == -1) || ...
        (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(m,jd)==-1) )        
	fprintf(fid, '"%s_std",', ar.model(m).data(jd).y);
end
fprintf(fid, '\n');

for jd = ar.model(m).plot(jplot).dLink
    fprintf('%i ', jd);
    for j=1:length(ar.model(m).data(jd).tExp)
        fprintnumtab(fid, ar.model(m).data(jd).tExp(j));
        
        % conditions
        if(~isempty(ar.model(m).data(jd).condition))
            for jp = 1:length(ar.model(m).data(jd).condition)
                fprintnumtab(fid, str2double(ar.model(m).data(jd).condition(jp).value));
            end
        end
        
        % y data
        for jj=1:size(ar.model(m).data(jd).yExp,2)
            if(ar.model(m).data(jd).logfitting(jj))
                fprintnumtab(fid, 10^ar.model(m).data(jd).yExp(j,jj));
            else
                fprintnumtab(fid, ar.model(m).data(jd).yExp(j,jj));
            end
        end
        
        % ystd data
        if( (ar.config.useFitErrorMatrix == 0 && ar.config.fiterrors == -1) || ...
                (ar.config.useFitErrorMatrix==1 && ar.config.fiterrors_matrix(m,jd)==-1) )  
            for jj=1:size(ar.model(m).data(jd).yExp,2)
                fprintnumtab(fid, ar.model(m).data(jd).yExpStd(j,jj));
            end
        end
        fprintf(fid, '\n');
    end
end

fclose(fid);
fprintf('done\n');

function fprintnumtab(fid, num, template)
if(~exist('template','var'))
    template = '"%s",';
end
strtmp = sprintf('%f', num);
strtmp = strrep(strtmp, ',', '.');
fprintf(fid, template, strtrim(strtmp));


