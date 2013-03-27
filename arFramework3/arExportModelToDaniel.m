function arExportModelToDaniel(m,d)

global ar

fid = fopen(sprintf('daniel_%s.csv',ar.model(m).name),'w');

fprintf(fid, '"Description","Rate"');
for jx = 1:length(ar.model(m).x)
    fprintf(fid, ',"%s"',ar.model(m).x{jx});
end
fprintf(fid, '\n');

for jv = 1:length(ar.model(m).fv)
    fprintf(fid, '"Reaktion%i","%s"', jv, ar.model(m).fv{jv});
    for jx = 1:length(ar.model(m).x)
        if(ar.model(m).N(jx,jv)~=0)
            fprintf(fid,',"%i"', ar.model(m).N(jx,jv));
        else
            fprintf(fid,',""');
        end
    end
    fprintf(fid, '\n');
end

fclose(fid);

fid = fopen(sprintf('daniel_%s_%s.csv',ar.model(m).name,ar.model(m).data(d).name),'w');

for jp = 1:length(ar.model(m).data(d).pold)
    fprintf(fid, '"%s","%s"\n', ar.model(m).data(d).pold{jp}, ar.model(m).data(d).fp{jp});
end

fclose(fid);