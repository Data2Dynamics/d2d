%  arImportSBML(filename, tEnd)
% 
%       tEnd        Default: 100
% 
% import SBML model and translate to .def files
% 
%   States and parameters consisting of a single character are replaced by
%   longer symbols.
%   State- and parameter names which coincide with mathematical functions
%   in symbolic the Symbolic Toolbox are replaced.
%   
% Example: 
% arImportSBML('BIOMD0000000379')
% 
% Example: 
%  ms = arImportSBML('BIOMD0000000379',100)
%   

function varargout = arImportSBML(filename, tEnd)
if(~exist('tEnd','var') || isempty(tEnd))
    tEnd = 100;
end

if exist('TranslateSBML','file')~=3
    warning('TranslateSBML not found. Please install libSBML and/or add it to Matlab''s search path. EXIT arImportSBML.m now.')
    return
end

try
    m = TranslateSBML([filename '.sbml']);
catch %#ok<CTCH>
    m = TranslateSBML([filename '.xml']);
end

mIn = m;

m = AdaptVariableNames(m);
m = rules2input(m);
m = findRateRules(m);

% [~,ia] = setdiff({m.species.id},{m.rule.variable});
% m.species = m.species(ia);  % Remove species which are combinations of real dynamic variables, in d2d something like observables

for i=1:length(m.compartment)
    if(sum(m.compartment(i).size<1e-5)>0)
        warning('Small compartment size exist. Could be rescale the equations because of integration problems.')
    end
end

%% model file
new_filename = strrep(filename,' ','_');
fid = fopen([new_filename '.def'], 'w');

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
    fprintf(fid, '%s\t T\t "%s"\t time\t 0\t %i\t\n', m.time_symbol, 'n/a', tEnd);
else
    fprintf(fid, 't\t T\t "%s"\t time\t 0\t %i\t\n', 'n/a', tEnd);
end

fprintf(fid, '\nCOMPARTMENTS\n');
for j=1:length(m.compartment)
    if(m.compartment(j).isSetSize)
        fprintf(fid, '%s\t V\t "%s"\t vol.\t %g\n', sym_check(m.compartment(j).id), 'n/a', ...
            m.compartment(j).size);
    else
        fprintf(fid, '%s\t V\t "%s"\t vol.\n', sym_check(m.compartment(j).id), 'n/a');
    end
end

fprintf(fid, '\nSTATES\n');
% for j=1:length(m.species)
pat = cell(0); % if length of species names ==1, the species names are extended by '_state'
rep = cell(0);
for j = find([m.species.isSetInitialAmount] | [m.species.isSetInitialConcentration] & ~[m.species.boundaryCondition]) % rules should not be defined as states, e.g. K_PP_norm in Huang1996 BIOMD0000000009
    if length(m.species(j).id)==1 %|| strcmp(m.species(j).id,'beta')==1  % special cases or too short
        pat{end+1} =  m.species(j).id; %#ok<AGROW>
        rep{end+1} = [m.species(j).id,'_state']; %#ok<AGROW>
        m.species(j).id2 = rep{end};        
        fprintf(fid, '%s\t C\t "%s"\t conc.\t %s\t 1\t "%s"\n', sym_check(rep{end}), 'n/a', ...
            m.species(j).compartment, m.species(j).name);
    else  % standard case
        m.species(j).id2 = m.species(j).id;
        fprintf(fid, '%s\t C\t "%s"\t conc.\t %s\t 1\t "%s"\n', sym_check(m.species(j).id), 'n/a', ...
            m.species(j).compartment, m.species(j).name);        
    end
end

for j=1:length(m.parameter)
    if length(m.parameter(j).id)==1 %|| strcmp(m.species(j).id,'beta')==1  % special cases or too short
        pat{end+1} =  m.parameter(j).id; %#ok<AGROW>
        rep{end+1} = [m.parameter(j).id,'_parameter']; %#ok<AGROW>
        m.parameter(j).id = rep{end};
    end
end  

for i=1:length(pat)
    for j=1:length(m.rule)
        m.rule(j).variable = mysubs(m.rule(j).variable,pat{i},rep{i});
        try
            m.rule(j).formula   = mysubs(m.rule(j).formula,pat{i},rep{i});
        catch
            m.rule(j).formula = regexprep(m.rule(j).formula,['(^|(\W)',pat{i},'($|\W)'],['$1',rep{i},'$2'],'all');
        end
    end
end

if(isfield(m,'raterule'))
    for i=1:length(pat)
        for j=1:length(m.raterule)
            m.raterule(j).variable = mysubs(m.raterule(j).variable,pat{i},rep{i});
            try
                m.raterule(j).formula   = mysubs(m.raterule(j).formula,pat{i},rep{i});
            catch
                m.raterule(j).formula = regexprep(m.raterule(j).formula,['(^|(\W)',pat{i},'($|\W)'],['$1',rep{i},'$2'],'all');
            end
        end
    end
end

fprintf(fid, '\nINPUTS\n');
for j=1:length(m.u)
    fprintf(fid, '%s\t C\t "%s"\t conc.\t%s\n', sym_check(m.u(j).variable), m.u(j).units, pow2mcode(m.u(j).formula,'power'));
end
% treat boundary species as constant inputs
for j=find([m.species.boundaryCondition])
    m.species(j).id2 = m.species(j).id;
    fprintf(fid, '%s\t C\t "%s"\t conc.\t"%s"\n', sym_check(m.species(j).id), 'n/a', ['init_' m.species(j).id]);
end

arWaitbar(0);
fprintf(fid, '\nREACTIONS\n');
if isfield(m,'raterule')
    for j=1:length(m.raterule)
        arWaitbar(j,length(m.raterule));
       
        prod_spec_name = m.raterule(j).variable;
        for i=1:length(rep)
            prod_spec_name = mysubs(prod_spec_name,pat{i},rep{i});
        end
        prod_spec_name = char(prod_spec_name);
        
        fprintf(fid,'\t -> %s', sym_check(prod_spec_name));

        tmpstr = sym(m.raterule(j).formula);
        % repace species names if too short
        for i=1:length(rep)
            tmpstr = mysubs(tmpstr,pat{i},rep{i});
        end
        
        tmpstr = replacePowerFunction(tmpstr);
        tmpstr = replacePowerFunction(tmpstr, 'pow');
        
        % replace rules
        tmpstr = sym(tmpstr);
        findrule = true;
        count = 0;
        while(findrule  && count<100)
            count = count+1;
            
            for jj=1:length(m.rule)
                try
                    tmpstr = mysubs(tmpstr, m.rule(jj).variable, ['(' m.rule(jj).formula ')']);
                catch
                    % rethrow(lasterr)
                    m.rule(jj).variable
                    m.rule(jj).formula
                end
                
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
        try
            tmpstr = replacePowerFunction(tmpstr);
            tmpstr = replacePowerFunction(tmpstr, 'pow');
        catch
            %         try
            tmpstr = pow2mcode(tmpstr);
            tmpstr = pow2mcode(tmpstr,'pow');
            %         end
        end
        
        fprintf(fid, ' \t CUSTOM "%s" \t"%s"\n', sym_check(tmpstr), m.raterule(j).name);

    end
else  % specified via reactions (standard case)
    for j=1:length(m.reaction)
        arWaitbar(j,length(m.reaction));
        for jj=1:length(m.reaction(j).reactant)
            % check if reactant is boundary species
            reactant_id = strcmp(m.reaction(j).reactant(jj).species,{m.species.id});
            isboundary = m.species(reactant_id).boundaryCondition;
            if ~isboundary
            react_spec_name = sym(m.reaction(j).reactant(jj).species);
                for i=1:length(rep)
                    react_spec_name = mysubs(react_spec_name,pat{i},rep{i});
                end
                react_spec_name = char(react_spec_name);

                if(~isnan(m.reaction(j).reactant(jj).stoichiometry))
                    stoichiometry = m.reaction(j).reactant(jj).stoichiometry;
                else
                    stoichiometry = 1;
                end
                for jjj=1:stoichiometry
                    fprintf(fid, '%s', sym_check(react_spec_name));
                    if(jj ~= length(m.reaction(j).reactant) || jjj ~= stoichiometry)
                        fprintf(fid, ' + ');
                    end
                end
            end
        end
        fprintf(fid, ' \t-> ');
        for jj=1:length(m.reaction(j).product)
            % check if product is boundary species
            product_id = strcmp(m.reaction(j).product(jj).species,{m.species.id});
            isboundary = m.species(product_id).boundaryCondition;
            if ~isboundary
                prod_spec_name = sym(m.reaction(j).product(jj).species);
                for i=1:length(rep)
                    prod_spec_name = mysubs(prod_spec_name,pat{i},rep{i});
                end
                prod_spec_name = char(prod_spec_name);

                if(~isnan(m.reaction(j).product(jj).stoichiometry))
                    stoichiometry = m.reaction(j).product(jj).stoichiometry;
                else
                    stoichiometry = 1;
                end
                for jjj=1:stoichiometry
                    fprintf(fid, '%s', sym_check(prod_spec_name));
                    if(jj ~= length(m.reaction(j).product) || jjj ~= stoichiometry)
                        fprintf(fid, ' + ');
                    end
                end
            end
        end
        
        tmpstr = sym(m.reaction(j).kineticLaw.math);
        % repace species names if too short
        for i=1:length(rep)
            tmpstr = mysubs(tmpstr,pat{i},rep{i});
        end
        
        % make parameters unique
        if(isfield(m.reaction(j).kineticLaw, 'parameter'))
            for jj=1:length(m.reaction(j).kineticLaw.parameter)
                tmpstr = mysubs(tmpstr, m.reaction(j).kineticLaw.parameter(jj).id, ...
                    [m.reaction(j).id '_' m.reaction(j).kineticLaw.parameter(jj).id]);
            end
        end
                
        % divide rates by compartment volume
        reaction_comp = findReactionCompartment(m,j);
        
        if ~isempty(reaction_comp) && sum(strcmp(reaction_comp,strsplit(char(tmpstr),'*')))==1
            tmpstr = [char(tmpstr) '/' sym_check(reaction_comp)];
            tmpstr = sym(tmpstr);
        end
            

        % % replace compartement volumes
        % for jj=1:length(m.compartment)
        %     tmpstr = mysubs(tmpstr, m.compartment(jj).id, num2str(m.compartment(jj).size));
        %     % tmpstr = mysubs(tmpstr, m.compartment(jj).id, 'vol_para');
        % end

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
        count = 0;
        while(findrule && count < 100)
            count = count+1;
            for jj=1:length(m.rule)
                tmpstr = mysubs(tmpstr, m.rule(jj).variable, ['(' m.rule(jj).formula ')']);
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
        try
            tmpstr = replacePowerFunction(tmpstr);
            tmpstr = replacePowerFunction(tmpstr, 'pow');
        catch
            %         try
            tmpstr = pow2mcode(tmpstr);
            tmpstr = pow2mcode(tmpstr,'pow');
            %         end
        end
     
        fprintf(fid, ' \t CUSTOM "%s" \t"%s"\n', sym_check(tmpstr), m.reaction(j).name);
    end
end   % end if dynamics specified either raterule or reaction

arWaitbar(-1);

fprintf(fid, '\nDERIVED\n');

fprintf(fid, '\nCONDITIONS\n');

fprintf(fid, '\nPARAMETERS\n');
for j=1:length(m.species)
    if(m.species(j).isSetInitialConcentration)
        ub = 1000;
        if(m.species(j).initialConcentration>ub)
            ub = m.species(j).initialConcentration*10;
        end
        fprintf(fid, 'init_%s\t %g\t %i\t 0\t 0\t %g\n', sym_check(m.species(j).id2), ...
            m.species(j).initialConcentration, m.species(j).constant==0, ub);
    elseif(m.species(j).isSetInitialAmount)
        ub = 1000;
        if(m.species(j).initialAmount>ub)
            ub = m.species(j).initialAmount*10;
        end
        fprintf(fid, 'init_%s\t %g\t %i\t 0\t 0\t %g\n', sym_check(m.species(j).id2), ...
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
            fprintf(fid, '%s\t %g\t %i\t 0\t 0\t %g\n', sym_check(m.parameter(j).id), ...
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

fid = fopen([new_filename '_data.def'], 'w');

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

if isfield(m.species,'id2')
    fprintf(fid, '\nOBSERVABLES\n');
    for j=1:length(m.species)    
        fprintf(fid, '%s_obs\t C\t "%s"\t conc.\t 0 0 "%s" "%s"\n', sym_check(m.species(j).id2), 'n/a', ...
            m.species(j).id2, m.species(j).name);
    end

    fprintf(fid, '\nERRORS\n');
    for j=1:length(m.species)
        fprintf(fid, '%s_obs\t "sd_%s"\n', sym_check(m.species(j).id2), sym_check(m.species(j).id2));
    end
end

fprintf(fid, '\nCONDITIONS\n');

fprintf(fid, '\nPARAMETERS\n');


fclose(fid);

if ~isdir('./Models')
    mkdir('Models');
end
system(['mv ',new_filename '.def Models']);

if ~isdir('./Data')
    mkdir('Data');
end
system(['mv ',new_filename '_data.def Data']);

if nargout > 0
    ms.d2d = m;
    ms.sbml = mIn;
    ms.pat = pat;
    ms.rep = rep;
    varargout{1} = ms;
end


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
    
    % Replace longest names first               %#<JV>
    [~,I]=sort(cellfun(@length,C), 'descend');  %#<JV>
    for j=1:length(D)
        pattern = sprintf('(^%s|(?<=[\\(\\+\\*\\-\\/])(%s)(?=[\\)\\+\\*\\-\\/]))', C{I(j)}, C{I(j)}); % use regex for replacing pars and vars correctly
        funtmplate = regexprep(funtmplate, pattern, ['(' D{I(j)} ')']);
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
    
    % Replace longest names first               %#<JV>
    [~,I]=sort(cellfun(@length,C), 'descend');  %#<JV>
    for j=1:length(D)
        funtmplate = strrep(funtmplate, C{I(j)}, ['(' D{I(j)} ')']); %#<JV> {j}=>I(j)
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


% Replaces power( ... , ...) syntax to matlab syntax (...)^(...)
% 
%   Examples:
% 
% powstr = 'k1 + power(power(k1*2, k2+(7*log(k3))),2) + 10*p3 + power(k1*2, k2+(7*log(k3))) + 10*p3';
% outstr = pow2mcode(powstr)
% 
% powstr2 = strrep(powstr,'power','pow');
% outstr = pow2mcode(powstr2,'pow')

function outstr = pow2mcode(powstr,funname)
if(~exist('funname','var') || isempty(funname))
    funname = 'power';
end

powstr = strrep(powstr,' ','');
outstr = '';

count = 0;
while strcmp(powstr,outstr)~=1 || count>50
    count = count+1;
    if(count>1)
        powstr = outstr;
    end
    outstr = pow2mcode_once(powstr,funname);
end

function outstr = pow2mcode_once(powstr,funname)

ind = regexp(powstr,[funname,'\('])+length(funname);

[kl,kr] = FindKlammerPaare(powstr,'()');
[kl,ia] = intersect(kl,ind);

if~isempty(ia)
    kr = kr(ia);
    [~,rf] = sort(kr-kl);
    
    tmp = powstr(kl(rf(1)):kr(rf(1)));
    ikomma = strfind(tmp,',');
    tmp = [tmp(1:ikomma-1),')^(',tmp(ikomma+1:end)];
    
    outstr = [powstr(1:kl(rf(1))-length(funname)-1),tmp,powstr(kr(rf(1))+1:end)];
else
    outstr = powstr;
end
    
    
% Findet in einer Formel (als string) die Paare "Klammer auf" kl und "Klammer zu" kr
% Bsp: ((...)(...))()
%  kl = 2     7     1    13
%  kr = 6    11    12    14
%  
% [kl,kr] = FindKlammerPaare(str)
% 
% [kl,kr] = FindKlammerPaare(str,klammern)
% 
%   klammer     String der L�nge 2
%     Default: klammern = '()'
%     oder z.B. '{}'

function [kl,kr] = FindKlammerPaare(str,klammern)
if(~exist('klammern','var') || isempty(klammern))
    klammern = '()';
end

indL = strfind(str,klammern(1));
indR = strfind(str,klammern(2));

if(isempty(indR))
	kr = [];
	kl = [];
else
    kr = NaN(size(indR));
    kl = NaN(size(indR));
	for i = 1:length(indR)
		kr(i) = indR(i);
		kl(i) = indL(max(find(indL<indR(i)))); %#ok<MXFND>
		indL(find(kl(i)==indL)) = []; %#ok<FNDSB>
	end
	% sort according to kl
	tmp = sortrows([kl;kr]')';
	kl = tmp(1,:);
	kr = tmp(2,:);
end

function m = findRateRules(m)

drin = [];
for i=1:length(m.rule)
    switch m.rule(i).typecode 
        case 'SBML_ASSIGNMENT_RULE'
            drin = [drin,i];%#ok<AGROW> % standard case
        case 'SBML_RATE_RULE'
            if ~isfield(m,'raterule')
                m.raterule = m.rule(i);
            else
                m.raterule(end+1) = m.rule(i);
            end
        otherwise
             m.rule(i).typecode
             error(' m.rule(i).typecode unknown');
    end

end
m.rule = m.rule(drin);


function m = rules2input(m)
is_input = zeros(size(m.rule));
for r=1:length(m.rule)
    s = symvar(m.rule(r).formula);
    if(~isempty(intersect(s,'TIME')))
        is_input(r) = 1;
    end
end

m.u = m.rule(find(is_input==1)); %#ok<FNDSB>
m.rule = m.rule(find(is_input~=1)); %#ok<FNDSB>

for i=1:length(m.u)
%     m.u(i).formula = char(mysubs(sym(m.u(i).formula),'time','t')); % does
%     not work, at least in R2014a
    m.u(i).formula = strrep(m.u(i).formula,'TIME','t');
end


% The following function will
%   alter variable names which cannot be handled by the Symbolic Math
%   function subs()

function m = AdaptVariableNames(m)

for i=1:length(m.species)
    m.species(i).id = sym_check(m.species(i).id);   
    m.species(i).compartment = sym_check(m.species(i).compartment);   
end


for i=1:length(m.parameter)
    m.parameter(i).id = sym_check(m.parameter(i).id);   
end

for i=1:length(m.rule)
    m.rule(i).variable = sym_check(m.rule(i).variable);   
    m.rule(i).formula = sym_check(m.rule(i).formula);   
end

for i=1:length(m.reaction)
    for j=1:length(m.reaction(i).reactant)
        m.reaction(i).reactant(j).species = sym_check(m.reaction(i).reactant(j).species);
    end
    m.reaction(i).kineticLaw.math = sym_check( m.reaction(i).kineticLaw.math); 
    if isfield(m.reaction(i).kineticLaw,'formula')
        m.reaction(i).kineticLaw.formula = sym_check( m.reaction(i).kineticLaw.formula);     
    end
end


function s = sym_check(s)
%Replacement if symbolic variable coincides with function. Without
%replacement subs would not work.
keywords = {'time','gamma','sin','cos','tan','beta','log','asin','atan','acos','acot','cot','theta','D','I','E'};

issym = strcmp(class(s),'sym'); %#ok<STISA>
if(issym)
    s = char(s);
end

sv = symvar(s);

svinter = intersect(sv,keywords);

for i=1:length(svinter)
    for j=1:3
        if(length(svinter{i})>1)
            s = regexprep(s,['(^|[\W])',svinter{i},'([\W]|$)'],['$1',upper(svinter{i}),'$2'],'all');
        else
            s = regexprep(s,['(^|[\W])',svinter{i},'([\W]|$)'],['$1',svinter{i},'_symbol','$2'],'all');
        end            
    end
end

if(issym)
    s = sym(s);
end


function s = mysubs(s,pat,rep)

keywords = {'time','gamma','sin','cos','tan','beta','log','asin','atan','acos','acot','cot','theta','D','I','E'};

issym = strcmp(class(s),'sym'); %#ok<STISA>
if(~issym)
   s = sym(s);
end
if isempty(intersect(pat,keywords))
    try
        s = char(subs(s,pat,rep));
        err=0;
    catch lasterr
        disp(lasterr)
        err=1;
    end
end
    
if ~isempty(intersect(pat,keywords)) || err==1
  % symbolic toolbox keywords (function) do not work with subs
    sv = symvar(char(s));
    svcell = cell(size(sv));
    for i=1:length(sv)
        svcell{i} = char(sv(i));
    end
    if ~isempty(intersect(pat,svcell))
        s = char(s);
        for j=1:3
            s = regexprep(s,['(^|(\W)',pat,'($|\W)'],['$1',rep,'$2'],'all');
        end
        if(issym)
            s = sym(s);
        end
    end
end


function c = findReactionCompartment(m, j)
% find compartment of reacting species to convert from SBML rate convention
% (particle flux) to d2d (concentration flux).

comp_r = {};
for jr=1:length(m.reaction(j).reactant);
    js = strcmp({m.species.id},m.reaction(j).reactant(jr).species);
    comp_r{jr} = m.species(js).compartment; %#ok<AGROW>
end

comp_p = {};
for jr=1:length(m.reaction(j).product);
    js = strcmp({m.species.id},m.reaction(j).product(jr).species);
    comp_p{jr} = m.species(js).compartment; %#ok<AGROW>
end

% check educt and product compartements for consistency
if ~isempty(comp_r)
    if length(unique(comp_r))~=1
        error('more than one educt compartment');
    end
end
if ~isempty(comp_p)
    if length(unique(comp_p))~=1
        error('more than one product compartment');
    end
end

% educt and product exist
if ~isempty(comp_r) && ~isempty(comp_p)
    % educt and product in the same compartment
    if isequal(unique(comp_r),unique(comp_p))
        c = comp_r{1};
    else
        c = [];
    end
% only educt has a compartment
elseif ~isempty(comp_r) && isempty(comp_p)
    c = comp_r{1};
% only product has a compartment
elseif isempty(comp_r) && ~isempty(comp_p)
    c = comp_p{1};
else
    c = [];
end

