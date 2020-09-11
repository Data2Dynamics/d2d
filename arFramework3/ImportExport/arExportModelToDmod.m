function filenames = arExportModelToDmod(whichone,m,d,prefix,replacements)
% arExportModelToDmod(whichone,m,d,prefix,replacements)
%
% Export model to dMod modeling framework
% (https://github.com/dkaschek/dMod) or ODESS analytical steady state
% solver (http://sysbio.uni-freiburg.de/mrosen/code9_release.py)
% 
%   whichone    Which kind of file, character or cell of characters
%               model       stoichometric matrix of the ODEs' righ-hand side
%               dataConditions the condition replacements specified in
%                           a data set (ar.model.data.pold and
%                           ar.model.data.fp) 
%               obs         observables of a data set (specified by d)
%               obsCalcUnits all dynamic variables and inputs as
%                           observables with same scaling parameter 
%                           for the same concentration unit
%               emptyobs    create an file which can be used as 2nd
%                           argument in the symmetry detection tool without
%                           the need to specify observations
% 
%   m           Model index (default: 1)
%   d           Data index  (default: 1:length(ar.model(m).data), not always required)
% 
%   prefix      Addionional prefix for the files.
%
%   replacements cell array of quantities that will be replaced with zeros
%   in reactions
% 
% Example (model and observables):
% arExportModelToDmod
% 
% Example (for unit calculation):
% arExportModelToDmod({'model','obsCalcUnits'})
% 
% Example:
% arExportModelToDmod('dataConditions')
% 
% Example:
% arExportModelToDmod('model')
% 
% Example:
% arExportModelToDmod({'model','emptyobs'})
%
% See also arExportDataToDaniel
global ar

if(~exist('whichone','var') || isempty(whichone))
    whichone = {'model','obs'};
end
if(~exist('prefix','var') || isempty(prefix))
    prefix = '';%sprintf('daniel_');
end
if(~exist('m','var') || isempty(m))
    m = 1;
end
if(~exist('d','var') || isempty(d))
    ds = 1:length(ar.model(m).data);
else
    ds = d;
end
if(~exist('replacements','var') || isempty(replacements))
    replacements = {''};
elseif ischar(replacements)
    replacements = {replacements};
end


if( iscell(whichone))
    filenames = cell(size(whichone));
    for i=1:length(whichone)
        filenames{i} = arExportModelToDmod(whichone{i},m,ds,prefix);
    end
    
else

    switch lower(whichone)
        case {'model'}
            filenames = sprintf('%s%s__model.csv',prefix,ar.model(m).name);
            fid = fopen(filenames,'w');
            
            fprintf(fid, '"Description","Rate"');
            for jx = 1:length(ar.model(m).x)
                fprintf(fid, ',"%s"',ar.model(m).x{jx});
            end
            fprintf(fid, '\n');
            
            for jv = 1:length(ar.model(m).fv)
                
                if isempty(replacements)
                    reaction = ar.model(m).fv{jv};
                else
                    reaction = ar.model(m).fv{jv};
                    for jr = 1:length(replacements)
                        reaction = char(arSubs(arSym(reaction), arSym(replacements{jr}), 0));
                    end
                end
                if reaction ~= '0'
                    fprintf(fid, '"Reaktion%i","%s"', jv, reaction);
                    for jx = 1:length(ar.model(m).x)
                        if(ar.model(m).N(jx,jv)~=0)
                            fprintf(fid,',"%i"', ar.model(m).N(jx,jv));
                        else
                            fprintf(fid,',""');
                        end
                    end
                fprintf(fid, '\n');
                end
            end            
            fclose(fid);            
        
        
        case {'data_conditions','dataconditions','dataconditions','datacondition','pold'} % alternative spelling
            for d=ds
                filenames = sprintf('%s%s__dataConditions__%s.csv',prefix,ar.model(m).name,ar.model(m).data(d).name);
                fid = fopen(filenames,'w');
                for jp = 1:length(ar.model(m).data(d).pold)
                    fprintf(fid, '"%s","%s"\n', ar.model(m).data(d).pold{jp}, ar.model(m).data(d).fp{jp});
                end
                fclose(fid);
            end
        
        case {'data','obs','observables'} % alternative spelling
            for d = ds
                d
                filenames = sprintf('%s%s__obs__%s.txt',prefix,ar.model(m).name,ar.model(m).data(d).name);
                fid = fopen(filenames,'w');
                
                for jd = 1:length(ar.model(m).data)
                    for jy = 1:length(ar.model(m).data(jd).y)
                        fprintf(fid, '%s',ar.model(m).data(jd).y{jy});
                        fprintf(fid, ' = %s \n',ar.model(m).data(jd).fy{jy});
                    end
                end
                fclose(fid);
            end
        
        case {'xobs','obscalcunits'} % alternative spelling
            filenames = sprintf('%s%s__obsCalcUnits.txt',prefix,ar.model(m).name);
            fid = fopen(filenames,'w');

            if(isempty(ar.model(m).uUnits))
                uuni = unique(ar.model(m).xUnits(:,1));
            elseif(isempty(ar.model(m).xUnits))
                uuni = unique(ar.model(m).uUnits(:,1));
            else
                uuni = unique([ar.model(m).xUnits(:,1);ar.model(m).uUnits(:,1)]);
            end
            for i=1:length(ar.model(m).x)
                unr = strmatch(ar.model(m).xUnits{i,1},uuni);
                fprintf(fid,'%s_obs = conc_unit%i*%s\n',ar.model(m).x{i},unr,ar.model(m).x{i});
            end
            for i=1:length(ar.model(m).u)
                unr = strmatch(ar.model(m).xUnits{i,1},uuni);
                fprintf(fid,'%s_obs = conc_unit%i*%s\n',ar.model(m).u{i},unr,ar.model(m).u{i});
            end
            fclose(fid);

        case {'emptyobs'}   % unfortunately the symmetry tool requires this argument, i.e. an empty file, if no observations are specified
            filenames = sprintf('%s%s__obs_empty.txt',prefix,ar.model(m).name);
            fid = fopen(filenames,'w');
            fprintf(fid,' ');
            fclose(fid);

        otherwise
            whichone
            error('Argument whichone unknown.');
    end
    

end
