function out = writeDynFunc(problem, funcfolder, modelname)
    %TODO
    
    funcname = strcat(modelname, '_dyn');
    filename = strcat(funcname, '.m');
    filepath = fullfile(funcfolder, filename);
    
    if isempty_ext(problem.amimodel)
        problem.getDynamics;
    end
    dx = string(problem.amimodel.xdot);
    
    fileId = fopen(filepath, 'w');
    
    fprintf(fileId, 'function dx = %s(t, x, p, c)\n', funcname);
    fprintf(fileId, '\tdx = zeros(%d, 1);\n\n', numel(dx));    
    
    for i = 1:numel(dx)
        dxi = parseExpression(problem, dx(i));
        fprintf(fileId, '\tdx(%d) = %s;\n', i, dxi);
    end
    
    fprintf(fileId, 'end');
    
    out = fclose(fileId);
end