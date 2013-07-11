% import SBML model and translate

function arImportSBML(filename)

try
    m = TranslateSBML([filename '.sbml']);
catch %#ok<CTCH>
    m = TranslateSBML([filename '.xml']);
end

%% model file

fid = fopen([filename '.def'], 'w');

fprintf(fid, 'DESCRIPTION\n');
if(~isempty(m.name))
    fprintf(fid, '"%s"', m.name);
    if(~isempty(m.id))
        fprintf(fid, '" (%s)"', m.id);
    end
    fprintf(fid, '\n');
else
    if(~isempty(m.id))
        fprintf(fid, '"%s"', m.id);
        fprintf(fid, '\n');
    end
end


fprintf(fid, '"SBML level %i version %i"\n', m.SBML_level, m.SBML_version);
if(~isempty(m.notes))
    fprintf(fid, '"%s"\n', m.notes);
end

fprintf(fid, '\nPREDICTOR\n');
if(~isempty(m.time_symbol))
    fprintf(fid, '%s\t T\t "%s"\t time\t 0\t %i\t\n', m.time_symbol, 'n/a', 100);
else
    fprintf(fid, 't\t T\t "%s"\t time\t 0\t %i\t\n', 'n/a', 100);
end

fprintf(fid, '\nCOMPARTMENTS\n');
for j=1:length(m.compartment)
    if(m.compartment(j).isSetSize)
        fprintf(fid, '%s\t V\t "%s"\t vol.\t %g\n', m.compartment(j).id, 'n/a', ...
            m.compartment(j).size);
    else
        fprintf(fid, '%s\t V\t "%s"\t vol.\n', m.compartment(j).id, 'n/a');
    end
end

fprintf(fid, '\nSTATES\n');
for j=1:length(m.species)
    fprintf(fid, '%s\t C\t "%s"\t conc.\t %s\t 1\t "%s"\n', m.species(j).id, 'n/a', ...
        m.species(j).compartment, m.species(j).name);
end

fprintf(fid, '\nINPUTS\n');
% for j=1:length(m.u)
%     fprintf(fid, '%s\t C\t "%s"\t conc.\t\n', m.u(j).ID, m.u(j).unit);
% end

arWaitbar(0);
fprintf(fid, '\nREACTIONS\n');
for j=1:length(m.reaction)
    arWaitbar(j,length(m.reaction));
    for jj=1:length(m.reaction(j).reactant)
        if(~isnan(m.reaction(j).reactant(jj).stoichiometry))
            stoichiometry = m.reaction(j).reactant(jj).stoichiometry;
        else
            stoichiometry = 1;
        end
        for jjj=1:stoichiometry
            fprintf(fid, '%s', m.reaction(j).reactant(jj).species);
            if(jj ~= length(m.reaction(j).reactant) || jjj ~= stoichiometry)
                fprintf(fid, ' + ');
            end
        end
    end
    fprintf(fid, ' \t-> ');
    for jj=1:length(m.reaction(j).product)
        if(~isnan(m.reaction(j).product(jj).stoichiometry))
            stoichiometry = m.reaction(j).product(jj).stoichiometry;
        else
            stoichiometry = 1;
        end
        for jjj=1:stoichiometry
            fprintf(fid, '%s', m.reaction(j).product(jj).species);
            if(jj ~= length(m.reaction(j).product) || jjj ~= stoichiometry)
                fprintf(fid, ' + ');
            end
        end
    end
    
    tmpstr = sym(m.reaction(j).kineticLaw.math);
    
    % make parameters unique
    if(isfield(m.reaction(j).kineticLaw, 'parameter'))
        for jj=1:length(m.reaction(j).kineticLaw.parameter)
            tmpstr = subs(tmpstr, m.reaction(j).kineticLaw.parameter(jj).id, ...
                [m.reaction(j).id '_' m.reaction(j).kineticLaw.parameter(jj).id]);
        end
    end
    
    % replace compartement volumes
    for jj=1:length(m.compartment)
        tmpstr = subs(tmpstr, m.compartment(jj).id, num2str(m.compartment(jj).size));
    end
    
    % remove functions
    tmpstr = char(tmpstr);
    for jj=1:length(m.functionDefinition)
        tmpfun = m.functionDefinition(jj).math;
        tmpfun = strrep(tmpfun, 'lambda(', '');
        tmpfun = tmpfun(1:end-1);
        
        C = textscan(tmpfun, '%s', 'Whitespace', ',');
        C = C{1};
        
        tmpstr = replaceFunction(tmpstr, m.functionDefinition(jj).id, C(1:end-1), C(end));
    end
    
    % replace power function
    tmpstr = replacePowerFunction(tmpstr);
    tmpstr = replacePowerFunction(tmpstr, 'pow');
   
    % replace rules
    tmpstr = sym(tmpstr);
    findrule = true;
    while(findrule)
        for jj=1:length(m.rule)
            tmpstr = subs(tmpstr, m.rule(jj).variable, ['(' m.rule(jj).formula ')']);
        end
        findrule = false;
        vars = symvar(tmpstr);
        for jj=1:length(m.rule)
            if(sum(ismember(vars, sym(m.rule(jj).variable)))>0) %R2013a compatible
                findrule = true;
            end
        end
    end
    tmpstr = char(tmpstr);
    
    % replace power function
    tmpstr = replacePowerFunction(tmpstr);
    tmpstr = replacePowerFunction(tmpstr, 'pow');
    
    fprintf(fid, ' \t CUSTOM "%s" \t"%s"\n', tmpstr, m.reaction(j).name);
end
arWaitbar(-1);

fprintf(fid, '\nINVARIANTS\n');

fprintf(fid, '\nCONDITIONS\n');

fprintf(fid, '\nPARAMETERS\n');
for j=1:length(m.species)
    if(m.species(j).isSetInitialConcentration)
        ub = 1000;
        if(m.species(j).initialConcentration>ub)
            ub = m.species(j).initialConcentration*10;
        end
        fprintf(fid, 'init_%s\t %g\t %i\t 0\t 0\t %g\n', m.species(j).id, ...
            m.species(j).initialConcentration, m.species(j).constant==0, ub);
    elseif(m.species(j).isSetInitialAmount)
        ub = 1000;
        if(m.species(j).initialAmount>ub)
            ub = m.species(j).initialAmount*10;
        end
        fprintf(fid, 'init_%s\t %g\t %i\t 0\t 0\t %g\n', m.species(j).id, ...
            m.species(j).initialAmount, m.species(j).constant==0, ub);
    end
end
for j=1:length(m.parameter)
    isRule = false;
    for jjj=1:length(m.rule)
        if(strcmp(m.rule(jjj).variable, m.parameter(j).id))
            isRule = true;
        end
        if(strcmp(m.rule(jjj).variable, m.parameter(j).name))
            isRule = true;
        end
    end
    if(~isRule)
        if(m.parameter(j).isSetValue)
            ub = 1000;
            if(m.parameter(j).value>ub)
                ub = m.parameter(j).value*10;
            end
            fprintf(fid, '%s\t %g\t %i\t 0\t 0\t %g\n', m.parameter(j).id, ...
                m.parameter(j).value, m.parameter(j).constant==0, ub);
        end
    end
end
% parameters from reactions
for j=1:length(m.reaction)
    if(isfield(m.reaction(j).kineticLaw, 'parameter'))
        for jj=1:length(m.reaction(j).kineticLaw.parameter)
            if(m.reaction(j).kineticLaw.parameter(jj).isSetValue)
                ub = 1000;
                if(m.reaction(j).kineticLaw.parameter(jj).value>ub)
                    ub = m.reaction(j).kineticLaw.parameter(jj).value*10;
                end
                fprintf(fid, '%s_%s\t %g\t %i\t 0\t 0\t %g\n', m.reaction(j).id, ...
                    m.reaction(j).kineticLaw.parameter(jj).id, ...
                    m.reaction(j).kineticLaw.parameter(jj).value, m.reaction(j).kineticLaw.parameter(jj).constant==0, ub);
            end
        end
    end
end

fclose(fid);

%% data file

fid = fopen([filename '_data.def'], 'w');

fprintf(fid, 'DESCRIPTION\n');
if(~isempty(m.name))
    fprintf(fid, '"data file for %s"\n', m.name);
end

fprintf(fid, '\nPREDICTOR\n');
if(~isempty(m.time_symbol))
    fprintf(fid, '%s\t T\t "%s"\t time\t 0\t %i\t\n', m.time_symbol, 'n/a', 100);
else
    fprintf(fid, 't\t T\t "%s"\t time\t 0\t %i\t\n', 'n/a', 100);
end

fprintf(fid, '\nINPUTS\n');

fprintf(fid, '\nOBSERVABLES\n');
for j=1:length(m.species)
    fprintf(fid, '%s_obs\t C\t "%s"\t conc.\t 0 0 "%s" "%s"\n', m.species(j).id, 'n/a', ...
        m.species(j).id, m.species(j).name);
end

fprintf(fid, '\nERRORS\n');
for j=1:length(m.species)
    fprintf(fid, '%s_obs\t "sd_%s"\n', m.species(j).id, m.species(j).id);
end

fprintf(fid, '\nINVARIANTS\n');

fprintf(fid, '\nCONDITIONS\n');

fprintf(fid, '\nRANDOM\n');

fprintf(fid, '\nPARAMETERS\n');


fclose(fid);



function str = replaceFunction(str, funstr, C, funmat)

% %% test replace functions
% clc
% str = 'k1 + power(k1*2, k2+(7*log(k3))) + 10*p3 + power(k1*2, k2+(7*log(k3))) + 10*p3';
% funstr = 'power';
% C{1} = 'a';
% C{2} = 'b';
% funmat = 'a^b';
% disp(str);

funindex = strfind(str, [funstr '(']);
while(~isempty(funindex))
    
    substr = str(funindex(1):end);
    
    openindex = strfind(substr, '(');
    closeindex = strfind(substr, ')');
    
    mergedindex = [openindex closeindex];
    rankingindex = [ones(size(openindex)) -ones(size(closeindex))];
    
    [sortedmergedindex, isortedindex] = sort(mergedindex);
    sortedrankingindex = rankingindex(isortedindex);
    
    endfunindex = find(cumsum(sortedrankingindex)==0);
    if(isempty(endfunindex))
        error('bracketing error close to function %s', funstr);
    end
    endfunindex = sortedmergedindex(endfunindex(1));
    
    substr = substr(openindex+1:endfunindex-1);
    
    D = textscan(substr, '%s', 'Whitespace', ',');
    D = D{1};
    if(length(C)~=length(D))
        error('input output parameter mismatch');
    end
    
    funtmplate = funmat;
    for j=1:length(D)
        funtmplate = strrep(funtmplate, C{j}, ['(' D{j} ')']);
    end
    funtmplate = ['(' funtmplate ')']; %#ok<AGROW>
    % disp(funtmplate)
    
    if(funindex(1)-1>1 && funindex(1)+endfunindex<length(str))
        str = [str(1:funindex(1)-1) funtmplate str(funindex(1)+endfunindex:end)];
    elseif(funindex(1)-1>1)
        str = [str(1:funindex(1)-1) funtmplate];  
    elseif(funindex(1)+endfunindex<length(str))
        str = [funtmplate str(funindex(1)+endfunindex:end)];
    else
        str = funtmplate;
    end
    str = cell2mat(str);
    % disp(str)
    
    funindex = strfind(str, funstr);
end
% disp(str)

str = char(sym(str));



function str = replacePowerFunction(str, funstr)

% %% test replace functions
% clc
% str = 'k1 + power(k1*2, k2+(7*log(k3))) + 10*p3 + power(k1*2, k2+(7*log(k3))) + 10*p3';

if(nargin<2)
    funstr = 'power';
end

C{1} = 'a';
C{2} = 'b';
funmat = 'a^b';

% disp(str);
funindex = strfind(str, [funstr '(']);
while(~isempty(funindex))
    
    substr = str(funindex(1):end);
    
    openindex = strfind(substr, '(');
    closeindex = strfind(substr, ')');
    
    mergedindex = [openindex closeindex];
    rankingindex = [ones(size(openindex)) -ones(size(closeindex))];
    
    [sortedmergedindex, isortedindex] = sort(mergedindex);
    sortedrankingindex = rankingindex(isortedindex);
    
    endfunindex = find(cumsum(sortedrankingindex)==0);
    if(isempty(endfunindex))
        error('bracketing error close to function %s', funstr);
    end
    endfunindex = sortedmergedindex(endfunindex(1));
    
    substr = substr(openindex+1:endfunindex-1);
    
    D = textscan(substr, '%s', 'Whitespace', ',');
    D = D{1};
    if(length(C)~=length(D))
        error('input output parameter mismatch');
    end
    
    funtmplate = funmat;
    for j=1:length(D)
        funtmplate = strrep(funtmplate, C{j}, ['(' D{j} ')']);
    end
    funtmplate = ['(' funtmplate ')']; %#ok<AGROW>
%     disp(funtmplate)
    
    if(funindex(1)-1>1 && funindex(1)-1+endfunindex<length(str)) % in between
        str = [str(1:funindex(1)-1) funtmplate str(funindex(1)+endfunindex:end)];
    elseif(funindex(1)-1>1) % at begining
        str = [str(1:funindex(1)-1) funtmplate];  
    elseif(funindex(1)-1+endfunindex<length(str)) % at end
        str = [funtmplate str(funindex(1)+endfunindex:end)];
    else % whole string
        str = funtmplate;
    end
%     disp(str)
    
    funindex = strfind(str, funstr);
end
% disp(str)

str = char(sym(str));

