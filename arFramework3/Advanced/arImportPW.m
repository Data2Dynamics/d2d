% arImportPW(filename)
%
% import PW model and translate
%
%   filename - name of def file

function arImportPW(filename)

m = eval(filename);

%% model file

fid = fopen([filename '.def'], 'w');

fprintf(fid, 'DESCRIPTION\n');
fprintf(fid, '"%s"\n', m.name); 
fprintf(fid, '"%s"\n', m.description);

fprintf(fid, '\nPREDICTOR\n');
fprintf(fid, 't\t T\t min\t time\t 0\t %i\t\n', max(m.t));

fprintf(fid, '\nCOMPARTMENTS\n');
for j=1:length(m.c)
    fprintf(fid, '%s\t V\t "%s"\t vol.\t %g\n', m.c(j).ID, m.c(j).unit, m.c(j).size);
end

fprintf(fid, '\nSTATES\n');
for j=1:length(m.x)
    fprintf(fid, '%s\t C\t "%s"\t conc.\t%s\n', m.x(j).ID, m.x(j).unit, m.x(j).compartment);
end

fprintf(fid, '\nINPUTS\n');
for j=1:length(m.u)
    fprintf(fid, '%s\t C\t "%s"\t conc.\t\n', m.u(j).ID, m.u(j).unit);
end

fprintf(fid, '\nREACTIONS\n');
for j=1:length(m.r)
    for jj=1:length(m.r(j).reactants)
        fprintf(fid, '%s', m.r(j).reactants{jj});
        if(jj < length(m.r(j).reactants))
            fprintf(fid, ' + ');
        else
            fprintf(fid, ' \t-> ');
        end
    end
    for jj=1:length(m.r(j).products)
        fprintf(fid, '%s', m.r(j).products{jj});
        if(jj < length(m.r(j).products))
            fprintf(fid, ' + ');
        end
    end
    
    tmpstr = m.r(j).rateSignature;
    for jj=length(m.r(j).reactants):-1:1
        tmpstr = strrep(tmpstr, sprintf('r%i', jj), m.r(j).reactants{jj});
    end
    for jj=length(m.r(j).products):-1:1
        tmpstr = strrep(tmpstr, sprintf('p%i', jj), m.r(j).products{jj});
    end
    for jj=length(m.r(j).modifiers):-1:1
        tmpstr = strrep(tmpstr, sprintf('m%i', jj), m.r(j).modifiers{jj});
    end
    for jj=length(m.r(j).parameters):-1:1
        tmpstr = strrep(tmpstr, sprintf('k%i', jj), m.r(j).parameters{jj});
    end
    
    fprintf(fid, ' \t CUSTOM "%s" \t"%s"\n', tmpstr, m.r(j).ID);
end

fprintf(fid, '\nINVARIANTS\n');

fprintf(fid, '\nCONDITIONS\n');
% for j=1:length(m.x)
%     if(strcmp(m.x(j).type,'fix'))
%         fprintf(fid, 'init_%s\t "%g"\n', m.x(j).ID, m.x(j).startValue);
%     end
% end
% for j=1:length(m.k)
%     if(strcmp(m.k(j).type,'fix'))
%         fprintf(fid, '%s\t "%g"\n', m.k(j).ID, m.k(j).value);
%     end
% end

fprintf(fid, '\nPARAMETERS\n');
for j=1:length(m.x)
    if(~strcmp(m.x(j).type,'fix'))
        fprintf(fid, 'init_%s\t %g\t 1\t 0\t %g\t %g\n', m.x(j).ID, m.x(j).startValue, ...
            m.x(j).minValue, m.x(j).maxValue);
    else
        fprintf(fid, 'init_%s\t %g\t 0\t 0\t %g\t %g\n', m.x(j).ID, m.x(j).startValue, ...
            m.x(j).minValue, m.x(j).maxValue);
    end
end
for j=1:length(m.k)
    if(~strcmp(m.k(j).type,'fix'))
        fprintf(fid, '%s\t %g\t 1\t 0\t %g\t %g\n', m.k(j).ID, m.k(j).value, ...
            m.k(j).minValue, m.k(j).maxValue);
    else
        fprintf(fid, '%s\t %g\t 0\t 0\t %g\t %g\n', m.k(j).ID, m.k(j).value, ...
            m.k(j).minValue, m.k(j).maxValue);
    end
end

fclose(fid);

%% data file

fid = fopen([filename '_data.def'], 'w');

fprintf(fid, 'DESCRIPTION\n');
fprintf(fid, '"%s"\n', m.name); 
fprintf(fid, '"%s"\n', m.description);

fprintf(fid, '\nPREDICTOR\n');
fprintf(fid, 't\t T\t min\t time\t 0\t %i\t\n', max(m.t));

fprintf(fid, '\nINPUTS\n');

fprintf(fid, '\nOBSERVABLES\n');
for j=1:length(m.y)
    tmpstr = m.y(j).scalingParameter;
    index = -1;
    for jj=1:length(m.s)
        if(strcmp(m.s(jj).ID, tmpstr))
            index = jj;
            break;
        end
    end
    if(index~=-1 && ~strcmp(m.s(index).type, 'fix'))
        fprintf(fid, '%s\t C\t "%s"\t conc.\t 0 0 "%s*(%s)"\n', m.y(j).ID, m.y(j).unit, tmpstr, m.y(j).rhsWithIDs);
    else
        fprintf(fid, '%s\t C\t "%s"\t conc.\t 0 0 "%s"\n', m.y(j).ID, m.y(j).unit, m.y(j).rhsWithIDs);
    end
end

fprintf(fid, '\nERRORS\n');
for j=1:length(m.y)
    fprintf(fid, '%s\t "sd_%s"\n', m.y(j).ID, m.y(j).ID);
%     if(strcmp(m.y(j).errorModel.type, 'rhs'))
%         fprintf(fid, '%s\t "%s"\n', m.y(j).ID, strrep(m.y(j).errorModel.rhs, 'y', m.y(j).ID));
%     else
%         error('error model not implemented');
%     end
end

fprintf(fid, '\nINVARIANTS\n');

fprintf(fid, '\nCONDITIONS\n');
% for j=1:length(m.s)
%     if(strcmp(m.s(j).type,'fix'))
%         fprintf(fid, '%s\t "%g"\n', m.s(j).ID, m.s(j).value);
%     end
% end

fprintf(fid, '\nRANDOM\n');

fprintf(fid, '\nPARAMETERS\n');
for j=1:length(m.s)
    if(strcmp(m.s(j).type,'fix'))
        fprintf(fid, '%s\t %g 0 0 \t%g \t%g\n', m.s(j).ID, m.s(j).value, ...
            m.s(j).minValue, m.s(j).maxValue);
    else
        fprintf(fid, '%s\t %g 1 0 \t%g \t%g\n', m.s(j).ID, m.s(j).value, ...
            m.s(j).minValue, m.s(j).maxValue);

    end
end

fclose(fid);
