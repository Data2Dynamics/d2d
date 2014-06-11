% create D2D file structure and template .def files

function arCreateProject

name = input('enter new project name: ', 's');

if(~isempty(name))
    path = ['./' name];
    
    if(exist(path, 'dir'))
        fprintf('folder %s already exists!\n', path);
    else
        mkdir(path)
        cd(path)
        
        mkdir('Models');
        cd('Models');
        copyfile(which('model_template.def'),'./model_template.def');
        cd ..
        
        mkdir('Data');
        cd('Data');
        copyfile(which('data_template.def'),'./data_template.def');
        copyfile(which('data_template.xls'),'./data_template.xls');
        cd ..
        
        fid = fopen('./Setup.m' , 'W');
        
        fprintf(fid, 'arInit;\n');
        fprintf(fid, 'arLoadModel(''model_template'');\n');
        fprintf(fid, 'arLoadData(''data_template'');\n');
        fprintf(fid, 'arCompileAll;\n');
        
        fclose(fid);
    end    
end