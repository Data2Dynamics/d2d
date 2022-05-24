function arStrikeGolddWriteModel


global ar

ppl = 6; % number of parameters per line
epl = 1; % number of equations per line

fprintf(' Writing the %s model ... ', ar.ia.modelname);

% start writing the model
fileID = fopen([ar.ia.modelname,'.m'],'w');
fprintf(fileID, '\nclear x h u p f ics known_ics\n\n');

% syms
fprintf(fileID, 'syms  ');
count=0;
for i=1:length(ar.ia.data.syms)
    fprintf(fileID, '%s ',ar.ia.data.syms{i});
    count = count + 1;
    if mod(count,ppl)==0
        fprintf(fileID, ' ...\n      ');
    end
end
fprintf(fileID, '\n\n\n');

% states
writeToFile('x', ar.ia.data.x, ppl, 'states', ';', fileID)

% outputs
writeToFile('h', ar.ia.data.h, epl, 'outputs', ';', fileID)

% known input
writeToFile('u', ar.ia.data.u, ppl, 'known inputs', ';', fileID)

% unknown input
writeToFile('w', ar.ia.data.w, ppl, 'unknown inputs', ';', fileID)

% parameters
writeToFile('p', ar.ia.data.p, ppl, 'parameters', ';', fileID)

% dynamic equations
writeToFile('f', ar.ia.data.f, epl, 'dynamic equations', ';', fileID)

% initial conditions
writeToFile('ics', ar.ia.data.ics, length(ar.ia.data.ics), 'initial conditions', ',', fileID)

% which initial conditions are known
writeToFile('known_ics', ar.ia.data.known_ics, length(ar.ia.data.known_ics), 'which initial conditions are known', ',', fileID)


% save
fprintf(fileID, '\n%% %s\n','save');
fprintf(fileID, 'save(''%s'',''x'',''h'',''u'',''w'',''p'',''f'',''ics'',''known_ics'');', ar.ia.modelname);
fprintf(fileID, '\n\n\n\n');

fclose(fileID);

fprintf('[done]');
fprintf('  (%s)\n', ar.ia.modelname_checkstr);

end



% write strings
function writeToFile(x, fx, ppl, des, sep, fileID)
fprintf(fileID, '\n%% %s\n',des);
fprintf(fileID, '%s = [',x);
count=0;
for i=1:length(fx)
    fprintf(fileID, '%s',fx{i});
    if i~=length(fx)
        fprintf(fileID, '%s ',sep);
    end
    count = count + 1;
    if mod(count,ppl)==0 && i~=length(fx)
        fprintf(fileID, '\n     ');
    end
end
fprintf(fileID, '];');
fprintf(fileID, '\n\n');
end
