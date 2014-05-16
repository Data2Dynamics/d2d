function arExportDatatoDaniel(m,d)

global ar;
fid = fopen(sprintf('daniel_%s_data.csv',ar.model(m).name),'w');

fprintf(fid, '"condition","name","time","value","uncertainty"');
fprintf(fid, '\n');
n_time=length(ar.model(m).data(d).tExp);

%for jx = 1:length(ar.model(m).condition) 
    for jN=1:length(ar.model(m).data(d).yNames)
        for jt=1:n_time
            fprintf(fid, '"dose1"');
            fprintf(fid, ',"%s"',ar.model(m).data(d).yNames{jN});
            fprintf(fid, ',"%i"',ar.model(m).data(d).tExp(jt));      
            if(ar.model(m).data(d).logfitting(jN)==1)
                
                fprintf(fid, ',"%i"',10^ar.model(m).data(d).yExp(jt,jN));
                if(ar.config.fiterrors==-1)
                    fprintf(fid, ',"%i"',10^ar.model(m).data(d).yExp(jt,jN)*ar.model(m).data(d).yExpstd(jt,jN));
                else
                    fprintf(fid, ',"%i"',10^ar.model(m).data(d).yExp(jt,jN)*ar.model(m).data(d).ystdExpSimu(jt,jN));
                end
            else
                
                fprintf(fid, ',"%i"',ar.model(m).data(d).yExp(jt,jN));
                if(ar.config.fiterrors==-1)
                    fprintf(fid, ',"%i"',ar.model(m).data(d).yExpstd(jt,jN));
                else
                    fprintf(fid, ',"%i"',ar.model(m).data(d).ystdExpSimu(jt,jN));
                end
                
            end
            fprintf(fid, '\n');
    
        end
    end 
%end

fclose(fid);

fid = fopen(sprintf('daniel_%s_info.csv',ar.model(m).name),'w');

fprintf(fid, '"condition info"');
fprintf(fid, '\n');
fprintf(fid, '\n');

for jC=1:length(ar.model(m).data(d).condition)
    fprintf(fid, '"%s"',ar.model(m).data(d).condition(jC).parameter);
    fprintf(fid, ',"%s"',ar.model(m).data(d).condition(jC).value);
    fprintf(fid, '\n');
end

fprintf(fid, '\n');
fprintf(fid, '"Time in %s"', ar.model(m).data(d).tUnits{2});
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '"Information about measurements"');
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '"Compartments","Scale"');
fprintf(fid, '\n');
for jC=1:length(ar.model(m).c)
    fprintf(fid, '"%s"',ar.model(m).c{jC});
    fprintf(fid, ',"%s"',ar.model(m).pc{jC});
    fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '"Name","Type","Units","log-scale","Fitted error","Compartment"');
fprintf(fid, '\n');
for jN=1:length(ar.model(m).data(d).yNames)
    fprintf(fid, '"%s"',ar.model(m).data(d).yNames{jN});
    fprintf(fid, ',"%s"',ar.model(m).data(d).yUnits{jN,3});
    fprintf(fid, ',"%s"',ar.model(m).data(d).yUnits{jN,2});
    fprintf(fid, ',"%i"',ar.model(m).data(d).logfitting(jN));
    fprintf(fid, ',"%i"',ar.config.fiterrors);
    fprintf(fid, ',"%s"',ar.model(m).c{ar.model(m).cLink(jN)});
    fprintf(fid, '\n');
end

fclose(fid);

end