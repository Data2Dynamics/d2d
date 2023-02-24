% arCreateProject()
% 
%    create D2D file structure and template .def files
% 
%    model_template   Which model def to choose, either one of the
%           following abbreviations or directly a model def name 
%           (see folder arFramework3/ProjectTemplate)
%           Default:  'model_template_ODEs.def'
% 
%           ode       the most basic ODE example 'model_template_ODEs.def'
%           ode_input ODEs with an input
%           abc       ODEs für A -> B -> C
%           hill      Hill functions
%           RTF       Retardet transient function

function arCreateProject(model_template)
if ~exist('model_template','var') || isempty(model_template)
    model_template = 'model_template_ODEs.def';
end

%% Handle choice of model and data template:
% replace if short name for model_template was used:
switch lower(model_template) 
    case {'ode','odes'} 
        model_template = 'model_template_ODEs.def';
    case {'ode_input','odes_input'} 
        model_template = 'model_template_ODEs_withInput.def';
    case {'abc'} 
        model_template = 'model_template_ABC.def';
    case {'hill','hillfunction'}
        model_template = 'model_template_HillFunctions.def';
    case {'tf'}
        model_template = 'model_template_TransientFunction.def';
    case {'tf2'}
        model_template = 'model_template_TransientFunction2.def';
    case {'rtf'}
        model_template = 'model_template_RetardedTransientFunction.def';
    otherwise
        % do nothing
end

% assignment from model template to data template:
setup_template = 'Setup_template_default.m'; % default, might be overwritten in the following:
switch model_template
    case 'model_template_ABC.def'
        data_template = 'data_template_ABC.def';
    case 'model_template_HillFunctions.def'
        data_template = 'data_template_HillFunctions.def';
    case 'model_template_TransientFunction.def'
        data_template = 'data_template_TransientFunction.def';   
        setup_template = 'Setup_template_TransientFunction.m';
    case 'model_template_TransientFunction2.def'
        data_template = 'data_template_TransientFunction2.def';   
        setup_template = 'Setup_template_TransientFunction.m';
    otherwise % choose the default template
        data_template = 'data_template.def';
end

fprintf('New project using %s as model template, %s as data template.\n',model_template,data_template);

%%
name = input('enter new project name: ', 's');
name = strrep(name,' ','_');

path0 = pwd;

try
    if(~isempty(name))
        path = ['./' name];
        
        if(exist(path, 'dir'))
            fprintf('folder %s already exists! Please choose another project name. \n', path);
        else
            mkdir(path)
            cd(path)
            
            mkdir('Models');
            cd('Models');
            w = which(model_template);
            if isempty(w)                
                error('model_template.def not found. Please check your path and/or execute arInit first.')
            end            
            copyfile(w,'./model_template.def');
            cd ..
            
            mkdir('Data');
            cd('Data');
            copyfile(which(data_template),'./data_template.def');
            copyfile(which(strrep(data_template,'.def','.xls')),'./data_template.xls');
            cd ..
            
            copyfile(which(setup_template),'./Setup.m');

            edit ./Setup.m
        end
    end
catch ERR
    cd(path0)
    try 
        fclose(fid)
    end
    rethrow(ERR)
end