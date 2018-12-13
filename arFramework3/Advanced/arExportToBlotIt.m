% arExportToBlotIt
%
% Exports current arStruct to csv to use for BlotIt in R
% i.e. Experiment;Condition;time;name;value
%
% Example:
%      arExportToBlotIt     (-> AllObservables_ForBlotIt.csv)
%      dat <- read.csv("AllObservables_ForBlotIt.csv",sep=";",dec=".")
%      result=alignME(data=dat, model="ys-sj",errmodel="sigmaR",fixed=ys~Condition, latent=sj~1, error=sigmaR~1,log=TRUE,normalize=T)
%      write.table(result,file="BlotIt_result.csv",quote=F,sep=";", row.names=FALSE,append=TRUE)
%
% https://github.com/dkaschek/blotIt2/blob/master/R/blotIt2.R

function arExportToBlotIt

global ar

%% find observables
obs = cell(0);
for m=1:length(ar.model)
    for d=1:length(ar.model(m).data)
        obs = unique(union(obs,ar.model(m).data(d).y));
    end
end

%% generate filenames
for o=1:length(obs)
    file{o} = ['ForBlotIt_',obs{o},'.csv'];
end

%% open files, write header
for o=1:length(obs)
    fid = fopen('ForBlotIt2_allObservables.csv','w');

    filehandle(o) = fopen(file{o},'w');
    fprintf('open %s ...\n',file{o});
%     fprintf(filehandle(o),'Index;Experiment;Condition;min;Target;Signal\n');    
    fprintf(filehandle(o),'Experiment;Condition;time;name;value\n');    
    
    fprintf(fid,'Experiment;Condition;time;name;value;ObsFun\n');        
end
anzrow = zeros(size(filehandle));


%% loop over all conditoins: collect data
tExp = cell(size(ar.model));
yExp = cell(size(ar.model));
yExpStd = cell(size(ar.model));
logscale = cell(size(ar.model));
fy = cell(size(ar.model));  % experiment index
for m=1:length(ar.model)
    tExp{m} = cell(size(ar.model(m).condition));
    yExp{m} = cell(size(ar.model(m).condition));
    yExpStd{m} = cell(size(ar.model(m).condition));
    fy{m} = cell(size(ar.model(m).condition));
    logscale{m} = cell(size(ar.model(m).condition));

    for c=1:length(ar.model(m).condition)
        tExp{m}{c} = cell(size(obs));
        yExp{m}{c} = cell(size(obs));
        yExpStd{m}{c} = cell(size(obs));
        fy{m}{c} = cell(size(obs));
        logscale{m}{c} = cell(size(obs));
        
        ds = ar.model(m).condition(c).dLink;
        for d=1:length(ds)
            y = ar.model(m).data(ds(d)).y;
            for iy=1:length(y)
                indo = strmatch(y{iy},obs,'exact'); % index within all observables
                if length(indo)~=1
                    error('Observable not found or found multiple times!? Programming error. Please check');
                end

                indnotnan = find(~isnan(ar.model(m).data(ds(d)).yExp(:,iy)));
                tExp{m}{c}{indo} = [tExp{m}{c}{indo}; ar.model(m).data(ds(d)).tExp(indnotnan)];
                logscale{m}{c}{indo} = ar.model(m).data(ds(d)).logfitting(iy);
                yExp{m}{c}{indo} = [yExp{m}{c}{indo}; ar.model(m).data(ds(d)).yExp(indnotnan,iy)];  % don't enforce non-log
                for it =1:length(indnotnan)
                    fy{m}{c}{indo} = [fy{m}{c}{indo}; ar.model(m).data(ds(d)).fy(iy) ];
                end
            end
        end
    end
end


row = 0;
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        for o=1:length(yExp{m}{c})
            anzt = sum(~isnan(yExp{m}{c}{o}));  % number of data points
            if anzt>5
                [~,~,fy_index] = unique(fy{m}{c}{o});
                for t=1:length(yExp{m}{c}{o})
                    if ~isnan(yExp{m}{c}{o}(t))
                        anzrow(o) = anzrow(o)+1;
                        row = row+1;
                        %                     fprintf(filehandle(o),'%i;%i;%i;%.2d;%d;%d\n',row,c,dExp{m}{c}{o}(t),tExp{m}{c}{o}(t),yExp{m}{c}{o}(t),yExpStd{m}{c}{o}(t));
%                         fprintf(filehandle(o),'%i;%i;%i;%.2d;%f\n',row,c,dExp{m}{c}{o}(t),tExp{m}{c}{o}(t),yExp{m}{c}{o}(t));
%                         fprintf(filehandle(o),'%i;%i;%i;%.2d;%s;%f\n',row,c,dExp{m}{c}{o}(t),tExp{m}{c}{o}(t),obs{o},yExp{m}{c}{o}(t));
                        %                     fprintf('%i;%i;%i;%.2d;%d;%d\n',row,c,dExp{m}{c}{o}(t),tExp{m}{c}{o}(t),yExp{m}{c}{o}(t),yExpStd{m}{c}{o}(t));
                        fprintf(filehandle(o),'%i;%i;%.2d;%s;%f;%s;%f\n',fy_index(t),c,tExp{m}{c}{o}(t),obs{o},yExp{m}{c}{o}(t),fy{m}{c}{o}{t},logscale{m}{c}{o});

                        fprintf(fid,'%i;%i;%.2d;%s;%f;%s;%f\n',fy_index(t),c,tExp{m}{c}{o}(t),obs{o},yExp{m}{c}{o}(t),fy{m}{c}{o}{t},logscale{m}{c}{o});

                    end
                end
            end
        end
    end
end



%% close files, write header
for o=1:length(obs)
    fclose(filehandle(o));
    if anzrow(o)==0
        system(['del ',file{o}]);
    end
end

fclose(fid);



