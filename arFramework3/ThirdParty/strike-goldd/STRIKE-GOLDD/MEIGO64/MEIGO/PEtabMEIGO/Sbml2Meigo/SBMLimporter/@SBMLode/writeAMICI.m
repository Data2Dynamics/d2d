function writeAMICI(obj, modelname, outdir, verbose)
    % writeAMICI writes the symbolic information from an SBMLode object
    % into an AMICI model definition file
    %
    % Parameters:
    %  modelname: target name of the model (_syms.m will be appended to the name )
    %
    % Return values:
    % void
    
    if nargin == 3
        verbose = false;
    end
    
    if verbose
        fprintf('writing file ...\n')
    end
    
    if ~exist(outdir)
        fprintf("Provided directory not found.")
        fid = fopen([modelname '_syms.m'],'w');
    else
        fid = fopen(fullfile(outdir, [modelname '_syms.m']),'w');
    end
    
    fprintf(fid,['function model = ' modelname '_syms()\n']);
    fprintf(fid,'\n');
    if(strcmp(obj.time_symbol,''))
        fprintf(fid,'t = sym(''t'');\n');
    else
        fprintf(fid,[obj.time_symbol ' = sym(''t'');\n']);
    end
    
    fprintf(fid,'\n');
    fprintf(fid,'avogadro = 6.02214179e23;');
    
%     fprintf(fid,'model.debug = true;\n');
    writeDefinition('STATES','x','state',obj,fid)
    writeDefinition('PARAMETERS','p','parameter',obj,fid)
    writeDefinition('CONDITIONS','k','condition',obj,fid)
    writeDerived('DYNAMICS','xdot','xdot',obj,fid)
    writeDerived('INITIALIZATION','x0','initState',obj,fid)
    writeDerived('OBSERVABLES','y','observable',obj,fid)
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,['% EVENTS\n']);
    for ievent = 1:length(obj.trigger)
        str_trigger = char(obj.trigger(ievent));
        str_bolus = strjoin(arrayfun(@char,obj.bolus(:,ievent),'UniformOutput',false),',');
        fprintf(fid,['model.event(' num2str(ievent) ') = amievent(' ...
            str_trigger ', ...\n' ...
            '[' str_bolus '], ...\n' ...
            '[]);\n']);
    end
    fprintf(fid,'\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    fprintf(fid,'function r = pow(x,y)\n');
    fprintf(fid,'\n');
    fprintf(fid,'    r = x^y;\n');
    fprintf(fid,'\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    fprintf(fid,'function r = power(x,y)\n');
    fprintf(fid,'\n');
    fprintf(fid,'    r = x^y;\n');
    fprintf(fid,'\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
    
    for ifun = 1:length(obj.funmath)
        fprintf(fid,['function r = ' obj.funarg{ifun} '\n']);
        fprintf(fid,'\n');
        fprintf(fid,['    r = ' obj.funmath{ifun} ';\n']);
        fprintf(fid,'\n');
        fprintf(fid,'end\n');
        fprintf(fid,'\n');
    end 
    
    for fun = {'factorial','cei','psi'}
        fprintUnsupportedFunctionError(fun{1},fid)
    end
    
    fclose(fid);
end

function writeDefinition(header,identifier,field,this,fid)
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,['%%%%\n%% ' header '\n']);
    if(length(this.(field))>0)
        vars = strjoin(cellfun(@char,num2cell(this.(field)),'UniformOutput',false));
        
        tmp_string = string(split(vars, ' '));
        tmp_idx = find(arrayfun(@exist, tmp_string));
        for i = 1:numel(tmp_idx)
            idx = tmp_idx(i);
            tmp = sprintf('%s = sym("%s");', tmp_string(idx), tmp_string(idx));
            fprintf(fid,[tmp '\n']);
        end
        
        tmp_string(tmp_idx) = [];
        vars = char(join(tmp_string, ' '));        
        
        fprintf(fid,['syms ' vars '\n']);
    end
    fprintf(fid,['model.' identifier ' = [' strjoin(cellfun(@char,num2cell(this.(field)),'UniformOutput',false),',') '];\n']);
end

function writeDerived(header,identifier,field,this,fid)
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,['%%%%\n%% ' header '\n']);
    fprintf(fid,'\n');
    if(strcmp(header,'OBSERVABLES'))
        fprintf(fid,['%% ' strjoin(cellfun(@char,num2cell(this.observable_name),'UniformOutput',false),'\n%% ')  '\n']);
    elseif (strcmp(header,'DYNAMICS'))
        compartment = string(this.compartment);
        compartment = char(join(compartment, ' '));
        fprintf(fid, ['syms ' compartment '\n']);
    end
    
    fprintf(fid,['model.' identifier ' = [' strjoin(cellfun(@char,num2cell(this.(field)),'UniformOutput',false),', ...\n') '];']);
end

function fprintUnsupportedFunctionError(functionName,fid)
    fprintf(fid,['function r = ' functionName '(x)\n']);
    fprintf(fid,'\n');
    fprintf(fid,['    error(''The ' functionName ' function is currently not supported!'');\n']);
    fprintf(fid,'\n');
    fprintf(fid,'end\n');
    fprintf(fid,'\n');
end