%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   fitLoaded                               %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Part of evaluate_CellLines, do fitting of observational parameters after
%  random receptor/surface levels are assigned
%
%
function fitLoaded(isRnd)
global ar

if(~exist('isRnd','var') || isempty(isRnd))
    isRnd = 1;
end

arLink
init_vec={'EGFR','ErbB2','ErbB3','IGF1R','Met'};
RTK_level=[];
for i=1:length(init_vec)
    RTK_level=[RTK_level strmatch(strcat('init_',init_vec{i}),ar.pLabel,'exact')];
end
inits=get_initals(ar.model(1).data(1).yExp(1,1:5),0);
if(isnan(inits))
    ar.model(1).data(1).yExp(1,1:5) = [rand*2.5 rand*2.5 (rand-0.3)*2.2 (rand-0.3)*2.2 (rand-0.3)*2.2];
    inits=get_initals(ar.model(1).data(1).yExp(1,1:5),0);
    if(isnan(inits))
        error('couldnt get good inits');
    end
end
ar.p(RTK_level)=inits;

ar.ub(RTK_level)=ar.p(RTK_level)+1;
ar.lb(RTK_level)=ar.p(RTK_level)-1;
ar.mean(RTK_level) = ar.p(RTK_level);
ar.type(RTK_level)=1;
ar.std(RTK_level)=0.05;

for i=6:length(ar.model(1).data(1).yNames)
    nr_tmp = strfind(ar.model(1).data(1).py(strncmp(ar.model(1).data(1).py,'offset_',7)),strrep(ar.model(1).data(1).yNames{i},'_au',''));
    nr_tmp = find(~cellfun(@isempty,nr_tmp));
    if(length(nr_tmp)>1)
        nr_tmp = nr_tmp(2);
    end
    nr_tmp = find(strcmp(ar.model(1).data(1).py(nr_tmp),ar.pLabel)==1);
    if(~isnan(ar.model(1).data(1).yExp(1,i)))
        ar.p(nr_tmp)=ar.model(1).data(1).yExp(1,i);
        if(ar.p(nr_tmp)>ar.ub(nr_tmp))
            ar.ub(nr_tmp)=ar.p(nr_tmp)+1;
        end
        if(ar.p(nr_tmp)<ar.lb(nr_tmp))
            ar.lb(nr_tmp)=ar.p(nr_tmp)-1;
        end   
    end
end
dFit = ar.model(1).plot(end).dLink;
for i=1:length(ar.model(1).data)
    if(sum(i==dFit)==0)
        ar.model(1).data(i).qFit(:) = 0;
        continue;
    end
    %empty=find(sum(isnan(ar.model(1).data(i).yExp))==length(ar.model(1).data(i).tExp));
    ar.model(1).data(i).qFit(:) = 1;
    ar.model(1).data(i).qFit(find(sum(isnan(ar.model(1).data(i).yExp))==length(ar.model(1).data(i).tExp)))=0;
    ar.model(1).data(i).qFit(~cellfun(@isempty,strfind(ar.model(1).data(i).yNames,'FACS_')))=0;
end

if(~isfield(ar,'qFit_bkp'))
    ar.qFit_bkp = ar.qFit;
end
ar.qFit(:) = 2;
matches = regexpi(ar.pLabel,'^(offset|scale)\w*H322M$');
if(~isRnd)
    ar.qFit(RTK_level) = 1;
end
ar.qFit(find(~cellfun(@isempty,matches))) = 1;
arFit

ar.model(1).qPlotYs(:)=0;
ar.model(1).qPlotYs(end)=1;