% Write function c-files
%
% arWriteCFiles(forcedCompile, debug_mode)
%
% forceFullCompile:   recompile all object files      [false]
% debug_mode:         debug mode for c code           [false]
%
% Copyright Andreas Raue 2011 (andreas.raue@fdm.uni-freiburg.de)

function arWriteCFiles(forcedCompile, debug_mode)

error(['The function arWriteCFiles is deprecated. Please use arCompileAll instead. ' ...
    'See https://github.com/Data2Dynamics/d2d/wiki/First%20steps for description of work flow.']);

warnreset = warning;
warning('off','symbolic:mupadmex:MuPADTextWarning');

global ar

if(isempty(ar))
	error('please initialize by arInit')
end

if(~exist('forcedCompile','var'))
	forcedCompile = false;
end
if(~exist('debug_mode','var'))
    debug_mode = false;
end

if(~exist([cd '/Compiled'], 'dir'))
	mkdir([cd '/Compiled'])
end
if(~exist([cd '/Compiled/' ar.info.c_version_code], 'dir'))
	mkdir([cd '/Compiled/' ar.info.c_version_code])
end

% Compiled folder hook for cluster usage
fid = fopen('./Compiled/arClusterCompiledHook.m', 'W');
fprintf(fid, 'function arClusterCompiledHook\n');
fclose(fid);

% Functions
fid = fopen(['./Compiled/' ar.info.c_version_code '/arSimuCalcFunctions.c'], 'W');

usePool = exist('gcp','file')>0 && ~isempty(gcp('nocreate'));

% model equations
fprintf('\n');
for m=1:length(ar.model)
    do_h = nan(1,length(ar.model(m).condition));
    do_c = nan(1,length(ar.model(m).condition));
    for c=1:length(ar.model(m).condition)
        if(~debug_mode)
            fprintf(fid, '#include "%s.h"\n', ar.model(m).condition(c).fkt);
        else
            fprintf(fid, '#include "%s.c"\n', ar.model(m).condition(c).fkt);
        end
        do_h(c) = forcedCompile || ~exist(['./Compiled/' ar.info.c_version_code '/' ar.model(m).condition(c).fkt '.h'],'file');
        do_c(c) = forcedCompile || ~exist(['./Compiled/' ar.info.c_version_code '/' ar.model(m).condition(c).fkt '.c'],'file');
    end
    
    c_version_code = ar.info.c_version_code;
    config = ar.config;
    model.name = ar.model(m).name;
    model.us = ar.model(m).us;
    model.xs = ar.model(m).xs;
    condition = ar.model(m).condition;
    if(usePool)
        parfor c=1:length(ar.model(m).condition)
            strtmp = '';
            if(do_h(c))
                fid_odeH = fopen(['./Compiled/' c_version_code '/' condition(c).fkt '.h'], 'W'); % create header file
                arWriteHFilesODE(fid_odeH, config, condition(c));
                fclose(fid_odeH);
                strtmp = 'header...';
            end
            if(do_c(c))
                fid_ode = fopen(['./Compiled/' c_version_code '/' condition(c).fkt '.c'], 'W');
                arWriteCFilesODE(fid_ode, config, model, condition(c));
                fclose(fid_ode);
                fprintf('writing condition m%i c%i, %s (%s)...%sdone\n', m, c, model.name, condition(c).checkstr, strtmp);
            else
                fprintf('writing condition m%i c%i, %s (%s)...%sskipped\n', m, c, model.name, condition(c).checkstr, strtmp);
            end
        end
    else
        for c=1:length(ar.model(m).condition)
            strtmp = '';
            if(do_h(c))
                fid_odeH = fopen(['./Compiled/' c_version_code '/' condition(c).fkt '.h'], 'W'); % create header file
                arWriteHFilesODE(fid_odeH, config, condition(c));
                fclose(fid_odeH);
                strtmp = 'header...';
            end
            if(do_c(c))
                fid_ode = fopen(['./Compiled/' c_version_code '/' condition(c).fkt '.c'], 'W');
                arWriteCFilesODE(fid_ode, config, model, condition(c));
                fclose(fid_ode);
                fprintf('writing condition m%i c%i, %s (%s)...%sdone\n', m, c, model.name, condition(c).checkstr, strtmp);
            else
                fprintf('writing condition m%i c%i, %s (%s)...%sskipped\n', m, c, model.name, condition(c).checkstr, strtmp);
            end
        end
    end
    
    if(isfield(ar.model(m), 'data'))
        do_h = nan(1,length(ar.model(m).data));
        do_d = nan(1,length(ar.model(m).data));
        for d=1:length(ar.model(m).data)
            if(~debug_mode)
                fprintf(fid, '#include "%s.h"\n', ar.model(m).data(d).fkt);
            else
                fprintf(fid, '#include "%s.c"\n', ar.model(m).data(d).fkt);
            end
            do_h(d) = forcedCompile || ~exist(['./Compiled/' ar.info.c_version_code '/' ar.model(m).data(d).fkt '.h'],'file');
            do_d(d) = forcedCompile || ~exist(['./Compiled/' ar.info.c_version_code '/' ar.model(m).data(d).fkt '.c'],'file');
        end
        
        data = ar.model(m).data;
        model_name = model.name;
        if(usePool)
            parfor d=1:length(ar.model(m).data)
                strtmp = '';
                if(do_h(d))
                    fid_obsH = fopen(['./Compiled/' c_version_code '/' data(d).fkt '.h'], 'W'); % create header file
                    arWriteHFilesOBS(fid_obsH, data(d));
                    fclose(fid_obsH);
                    strtmp = 'header...';
                end
                if(do_d(d))
                    fid_obs = fopen(['./Compiled/' c_version_code '/' data(d).fkt '.c'], 'W');
                    arWriteCFilesOBS(fid_obs, config, data(d));
                    fclose(fid_obs);
                    fprintf('writing data m%i d%i, %s (%s)...%sdone\n', m, d, model_name, data(d).checkstr, strtmp);
                else
                    fprintf('writing data m%i d%i, %s (%s)...%sskipped\n', m, d, model_name, data(d).checkstr, strtmp);
                end
            end
        else
            for d=1:length(ar.model(m).data)
                strtmp = '';
                if(do_h(d))
                    fid_obsH = fopen(['./Compiled/' c_version_code '/' data(d).fkt '.h'], 'W'); % create header file
                    arWriteHFilesOBS(fid_obsH, data(d));
                    fclose(fid_obsH);
                    strtmp = 'header...';
                end
                if(do_d(d))
                    fid_obs = fopen(['./Compiled/' c_version_code '/' data(d).fkt '.c'], 'W');
                    arWriteCFilesOBS(fid_obs, config, data(d));
                    fclose(fid_obs);
                    fprintf('writing data m%i d%i, %s (%s)...%sdone\n', m, d, model_name, data(d).checkstr, strtmp);
                else
                    fprintf('writing data m%i d%i, %s (%s)...%sskipped\n', m, d, model_name, data(d).checkstr, strtmp);
                end
            end
        end
    end
end
fprintf(fid, '\n');

% map CVodeInit to fx
fprintf(fid, ' int AR_CVodeInit(void *cvode_mem, N_Vector x, double t, int im, int ic){\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if(im==%i & ic==%i) return CVodeInit(cvode_mem, fx_%s, RCONST(t), x);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '  return(-1);\n');
fprintf(fid, '}\n\n');

% map fx
fprintf(fid, ' void fx(realtype t, N_Vector x, double *xdot, void *user_data, int im, int ic){\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if(im==%i & ic==%i) fxdouble_%s(t, x, xdot, user_data);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '}\n\n');

% map fx0
fprintf(fid, ' void fx0(N_Vector x0, void *user_data, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if(im==%i & ic==%i) fx0_%s(x0, data);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '}\n\n');

% map CVDlsSetDenseJacFn to dfxdx
fprintf(fid, ' int AR_CVDlsSetDenseJacFn(void *cvode_mem, int im, int ic){\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if(im==%i & ic==%i) return CVDlsSetDenseJacFn(cvode_mem, dfxdx_%s);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '  return(-1);\n');
fprintf(fid, '}\n\n');

% map fsx0
fprintf(fid, ' void fsx0(int is, N_Vector sx_is, void *user_data, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '  if(im==%i & ic==%i) fsx0_%s(is, sx_is, data);\n', ...
            m-1, c-1, ar.model(m).condition(c).fkt);
    end
end
fprintf(fid, '}\n\n');

% map CVodeSensInit1 to fsx
fprintf(fid, ' int AR_CVodeSensInit1(void *cvode_mem, int nps, int sensi_meth, int sensirhs, N_Vector *sx, int im, int ic){\n');
if(ar.config.useSensiRHS)
    fprintf(fid, '  if (sensirhs == 1) {\n');
    for m=1:length(ar.model)
        for c=1:length(ar.model(m).condition)
            fprintf(fid, '    if(im==%i & ic==%i) return CVodeSensInit1(cvode_mem, nps, sensi_meth, fsx_%s, sx);\n', ...
                m-1, c-1, ar.model(m).condition(c).fkt);
        end
    end
    fprintf(fid, '  } else {\n');
end
for m=1:length(ar.model)
    for c=1:length(ar.model(m).condition)
        fprintf(fid, '    if(im==%i & ic==%i) return CVodeSensInit1(cvode_mem, nps, sensi_meth, NULL, sx);\n', ...
            m-1, c-1);
    end
end
if(ar.config.useSensiRHS)
    fprintf(fid, '  }\n');
end
fprintf(fid, '  return(-1);\n');
fprintf(fid, '}\n\n');

% map fu
fprintf(fid, ' void fu(void *user_data, double t, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if(im==%i & ic==%i) fu_%s(data, t);\n', ...
			m-1, c-1, ar.model(m).condition(c).fkt);
	end
end
fprintf(fid, '}\n\n');

% map fsu
fprintf(fid, ' void fsu(void *user_data, double t, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if(im==%i & ic==%i) fsu_%s(data, t);\n', ...
			m-1, c-1, ar.model(m).condition(c).fkt);
	end
end
fprintf(fid, '}\n\n');

% map fv
fprintf(fid, ' void fv(void *user_data, double t, N_Vector x, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if(im==%i & ic==%i) fv_%s(t, x, data);\n', ...
			m-1, c-1, ar.model(m).condition(c).fkt);
	end
end
fprintf(fid, '}\n\n');

% map fsv
fprintf(fid, ' void fsv(void *user_data, double t, N_Vector x, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if(im==%i & ic==%i) {\n\tdvdp_%s(t, x, data);\n\tdvdu_%s(t, x, data);\n\tdvdx_%s(t, x, data);\n}\n', ...
			m-1, c-1, ar.model(m).condition(c).fkt, ar.model(m).condition(c).fkt, ar.model(m).condition(c).fkt);
	end
end
fprintf(fid, '}\n\n');

% map dfxdp
fprintf(fid, ' void dfxdp(void *user_data, double t, N_Vector x, double *dfxdp, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if(im==%i & ic==%i) dfxdp_%s(t, x, dfxdp, data);\n', ...
			m-1, c-1, ar.model(m).condition(c).fkt);
	end
end
fprintf(fid, '}\n\n');

% map fy
fprintf(fid, ' void fy(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int iruns, double *y, double *p, double *u, double *x, int im, int id){\n');
for m=1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            fprintf(fid, '  if(im==%i & id==%i) fy_%s(t, nt, it, ntlink, itlink, ny, nx, iruns, y, p, u, x);\n', ...
                m-1, d-1, ar.model(m).data(d).fkt);
        end
    end
end
fprintf(fid, '}\n\n');

% map fystd
fprintf(fid, ' void fystd(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x, int im, int id){\n');
for m=1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            fprintf(fid, '  if(im==%i & id==%i) fystd_%s(t, nt, it, ntlink, itlink, ystd, y, p, u, x);\n', ...
                m-1, d-1, ar.model(m).data(d).fkt);
        end
    end
end
fprintf(fid, '}\n\n');

% map fsy
fprintf(fid, ' void fsy(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *su, double *sx, int im, int id){\n');
for m=1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            fprintf(fid, '  if(im==%i & id==%i) fsy_%s(t, nt, it, ntlink, itlink, sy, p, u, x, su, sx);\n', ...
                m-1, d-1, ar.model(m).data(d).fkt);
        end
    end
end
fprintf(fid, '}\n\n');

% map fsystd
fprintf(fid, ' void fsystd(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *sy, double *su, double *sx, int im, int id){\n');
for m=1:length(ar.model)
    if(isfield(ar.model(m), 'data'))
        for d=1:length(ar.model(m).data)
            fprintf(fid, '  if(im==%i & id==%i) fsystd_%s(t, nt, it, ntlink, itlink, systd, p, y, u, x, sy, su, sx);\n', ...
                m-1, d-1, ar.model(m).data(d).fkt);
        end
    end
end
fprintf(fid, '}\n\n');

% for arSSACalc
% call to fu and fv
fprintf(fid, '/* for arSSACalc.c */\n\n');
fprintf(fid, ' void fvSSA(void *user_data, double t, N_Vector x, int im, int ic){\n');
fprintf(fid, '  UserData data = (UserData) user_data;\n');
for m=1:length(ar.model)
	for c=1:length(ar.model(m).condition)
		fprintf(fid, '  if(im==%i & ic==%i) {\n', m-1, c-1);
        fprintf(fid, '    fu_%s(data, t);\n', ar.model(m).condition(c).fkt);
        fprintf(fid, '    fv_%s(t, x, data);\n', ar.model(m).condition(c).fkt);
        fprintf(fid, '  }\n');
	end
end
fprintf(fid, '}\n\n');

fclose(fid);

% compile
arCompile(forcedCompile);

% refresh file cache
rehash

warning(warnreset);


% write ODE headers
function arWriteHFilesODE(fid, config, condition)

fprintf(fid, '#ifndef _MY_%s\n', condition.fkt);
fprintf(fid, '#define _MY_%s\n\n', condition.fkt);

fprintf(fid, '#include <cvodes/cvodes.h>\n'); 
fprintf(fid, '#include <cvodes/cvodes_dense.h>\n');
fprintf(fid, '#include <nvector/nvector_serial.h>\n');
fprintf(fid, '#include <sundials/sundials_types.h>\n'); 
fprintf(fid, '#include <sundials/sundials_math.h>\n');  
fprintf(fid, '#include <udata.h>\n');
fprintf(fid, '#include <math.h>\n');
fprintf(fid, '#include <mex.h>\n');
fprintf(fid, '#include <arInputFunctionsC.h>\n');
fprintf(fid,'\n\n\n');

fprintf(fid, ' void fu_%s(void *user_data, double t);\n', condition.fkt);
fprintf(fid, ' void fsu_%s(void *user_data, double t);\n', condition.fkt);
fprintf(fid, ' void fv_%s(realtype t, N_Vector x, void *user_data);\n', condition.fkt);
fprintf(fid, ' void dvdx_%s(realtype t, N_Vector x, void *user_data);\n', condition.fkt);
fprintf(fid, ' void dvdu_%s(realtype t, N_Vector x, void *user_data);\n', condition.fkt);
fprintf(fid, ' void dvdp_%s(realtype t, N_Vector x, void *user_data);\n', condition.fkt);
fprintf(fid, ' int fx_%s(realtype t, N_Vector x, N_Vector xdot, void *user_data);\n', condition.fkt);
fprintf(fid, ' void fxdouble_%s(realtype t, N_Vector x, double *xdot_tmp, void *user_data);\n', condition.fkt);
fprintf(fid, ' void fx0_%s(N_Vector x0, void *user_data);\n', condition.fkt);
% fprintf(fid, ' int dfxdx_%s(int N, realtype t, N_Vector x,', condition.fkt); % sundials 2.4.0
fprintf(fid, ' int dfxdx_%s(long int N, realtype t, N_Vector x,', condition.fkt); % sundials 2.5.0
fprintf(fid, 'N_Vector fx, DlsMat J, void *user_data,');
fprintf(fid, 'N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);\n');
if(config.useSensiRHS)
    fprintf(fid, ' int fsx_%s(int Ns, realtype t, N_Vector x, N_Vector xdot,', condition.fkt);
    fprintf(fid, 'int ip, N_Vector sx, N_Vector sxdot, void *user_data,');
    fprintf(fid, 'N_Vector tmp1, N_Vector tmp2);\n');
end
fprintf(fid, ' void fsx0_%s(int ip, N_Vector sx0, void *user_data);\n', condition.fkt);
fprintf(fid, ' void dfxdp_%s(realtype t, N_Vector x, double *dfxdp, void *user_data);\n\n', condition.fkt);
fprintf(fid, '#endif /* _MY_%s */\n', condition.fkt);

fprintf(fid,'\n\n\n');


% ODE
function arWriteCFilesODE(fid, config, model, condition)

timedebug = false;

fprintf(fid, '#include "%s.h"\n',  condition.fkt);
fprintf(fid, '#include <cvodes/cvodes.h>\n');    
fprintf(fid, '#include <cvodes/cvodes_dense.h>\n');
fprintf(fid, '#include <nvector/nvector_serial.h>\n');
fprintf(fid, '#include <sundials/sundials_types.h>\n'); 
fprintf(fid, '#include <sundials/sundials_math.h>\n');  
fprintf(fid, '#include <udata.h>\n');
fprintf(fid, '#include <math.h>\n');
fprintf(fid, '#include <mex.h>\n');
fprintf(fid, '#include <arInputFunctionsC.h>\n');
fprintf(fid,'\n\n\n');

% write fu
fprintf(fid, ' void fu_%s(void *user_data, double t)\n{\n', condition.fkt);
if(timedebug) 
    fprintf(fid, '  printf("%%g \\t fu\\n", t);\n'); 
end;
if(~isempty(model.us))
    fprintf(fid, '  UserData data = (UserData) user_data;\n');
    fprintf(fid, '  double *p = data->p;\n');
    writeCcode(fid, condition, 'fu');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write fsu
fprintf(fid, ' void fsu_%s(void *user_data, double t)\n{\n', condition.fkt);
if(config.useSensis)
    if(sum(logical(condition.sym.dfudp(:)~=0))>0)
        fprintf(fid, '  UserData data = (UserData) user_data;\n');
        fprintf(fid, '  double *p = data->p;\n');
    
        writeCcode(fid, condition, 'fsu');
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write v
fprintf(fid, ' void fv_%s(realtype t, N_Vector x, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
    fprintf(fid, '  printf("%%g \\t fv\\n", t);\n');
end
if(~isempty(model.xs))
    fprintf(fid, '  UserData data = (UserData) user_data;\n');
    fprintf(fid, '  double *p = data->p;\n');
    fprintf(fid, '  double *u = data->u;\n');
    fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
    writeCcode(fid, condition, 'fv');
end

fprintf(fid, '\n  return;\n}\n\n\n');

% write dvdx
fprintf(fid, ' void dvdx_%s(realtype t, N_Vector x, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
    fprintf(fid, '  printf("%%g \\t dvdx\\n", t);\n');
end
if(~isempty(model.xs))
    fprintf(fid, '  UserData data = (UserData) user_data;\n');
    fprintf(fid, '  double *p = data->p;\n');
    fprintf(fid, '  double *u = data->u;\n');
    fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
    writeCcode(fid, condition, 'dvdx');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write dvdu
fprintf(fid, ' void dvdu_%s(realtype t, N_Vector x, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
    fprintf(fid, '  printf("%%g \\t dvdu\\n", t);\n');
end
if(~isempty(model.us) && ~isempty(model.xs))
    fprintf(fid, '  UserData data = (UserData) user_data;\n');
    fprintf(fid, '  double *p = data->p;\n');
    fprintf(fid, '  double *u = data->u;\n');
    fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
    writeCcode(fid, condition, 'dvdu');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write dvdp
fprintf(fid, ' void dvdp_%s(realtype t, N_Vector x, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
	fprintf(fid, '  printf("%%g \\t dvdp\\n", t);\n');
end
if(~isempty(model.xs))
    fprintf(fid, '  UserData data = (UserData) user_data;\n');
    fprintf(fid, '  double *p = data->p;\n');
    fprintf(fid, '  double *u = data->u;\n');
    fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
    if(~isempty(condition.sym.dfvdp))
        writeCcode(fid, condition, 'dvdp');
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write fx
fprintf(fid, ' int fx_%s(realtype t, N_Vector x, N_Vector xdot, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
    fprintf(fid, '  printf("%%g \\t fx\\n", t);\n');
end
if(~isempty(model.xs))
    fprintf(fid, '  int is;\n');
    fprintf(fid, '  UserData data = (UserData) user_data;\n');
    fprintf(fid, '  double *qpositivex = data->qpositivex;\n');
    fprintf(fid, '  double *p = data->p;\n');
    fprintf(fid, '  double *u = data->u;\n');
    fprintf(fid, '  double *v = data->v;\n');
    fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
    fprintf(fid, '  double *xdot_tmp = N_VGetArrayPointer(xdot);\n');
    fprintf(fid, '  fu_%s(data, t);\n', condition.fkt);
    fprintf(fid, '  fv_%s(t, x, data);\n', condition.fkt);
    writeCcode(fid, condition, 'fx');
    fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(model.xs));
    fprintf(fid, '    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;\n');
    fprintf(fid, '    if(qpositivex[is]>0.5 && x_tmp[is]<0.0 && xdot_tmp[is]<0.0) xdot_tmp[is] = -xdot_tmp[is];\n');
    fprintf(fid, '  }\n');
end

fprintf(fid, '\n  return(0);\n}\n\n\n');

% write fxdouble
fprintf(fid, ' void fxdouble_%s(realtype t, N_Vector x, double *xdot_tmp, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
    fprintf(fid, '  printf("%%g \\t fxdouble\\n", t);\n');
end
if(~isempty(model.xs))
    fprintf(fid, '  int is;\n');
    fprintf(fid, '  UserData data = (UserData) user_data;\n');
    fprintf(fid, '  double *p = data->p;\n');
    fprintf(fid, '  double *u = data->u;\n');
    fprintf(fid, '  double *v = data->v;\n');
    fprintf(fid, '  fu_%s(data, t);\n', condition.fkt);
    fprintf(fid, '  fv_%s(t, x, data);\n', condition.fkt);
    writeCcode(fid, condition, 'fx');
    fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(model.xs));
    fprintf(fid, '    if(mxIsNaN(xdot_tmp[is])) xdot_tmp[is] = 0.0;\n');
    fprintf(fid, '  }\n');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write fx0
fprintf(fid, ' void fx0_%s(N_Vector x0, void *user_data)\n{\n', condition.fkt);
if(~isempty(model.xs))
    fprintf(fid, '  UserData data = (UserData) user_data;\n');
    fprintf(fid, '  double *p = data->p;\n');
    fprintf(fid, '  double *u = data->u;\n');
    fprintf(fid, '  double *x0_tmp = N_VGetArrayPointer(x0);\n');
    writeCcode(fid, condition, 'fx0');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write dfxdx
% fprintf(fid, ' int dfxdx_%s(int N, realtype t, N_Vector x, \n', condition.fkt); % sundials 2.4.0
fprintf(fid, ' int dfxdx_%s(long int N, realtype t, N_Vector x, \n', condition.fkt); % sundials 2.5.0
fprintf(fid, '  \tN_Vector fx, DlsMat J, void *user_data, \n');
fprintf(fid, '  \tN_Vector tmp1, N_Vector tmp2, N_Vector tmp3)\n{\n');
if(timedebug)
    fprintf(fid, '  printf("%%g \\t dfxdx\\n", t);\n');
end

if(~isempty(model.xs))
    if(config.useJacobian)
        fprintf(fid, '  int is;\n');
        fprintf(fid, '  UserData data = (UserData) user_data;\n');
        fprintf(fid, '  double *p = data->p;\n');
        fprintf(fid, '  double *u = data->u;\n');
        fprintf(fid, '  double *dvdx = data->dvdx;\n');
        % fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
        fprintf(fid, '  dvdx_%s(t, x, data);\n', condition.fkt);
        writeCcode(fid, condition, 'dfxdx');
        fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(model.xs)^2);
        fprintf(fid, '    if(mxIsNaN(J->data[is])) J->data[is] = 0.0;\n');
        fprintf(fid, '  }\n');
    end
end
fprintf(fid, '\n  return(0);\n}\n\n\n');

% write fsv & fsx
if(config.useSensiRHS)
    fprintf(fid, ' int fsx_%s(int Ns, realtype t, N_Vector x, N_Vector xdot, \n', condition.fkt);
    fprintf(fid, '  \tint ip, N_Vector sx, N_Vector sxdot, void *user_data, \n');
    fprintf(fid, '  \tN_Vector tmp1, N_Vector tmp2)\n{\n');
    
    if(~isempty(model.xs))
        if(config.useSensis)
            fprintf(fid, '  int is;\n');
            fprintf(fid, '  UserData data = (UserData) user_data;\n');
            fprintf(fid, '  double *p = data->p;\n');
            fprintf(fid, '  double *u = data->u;\n');
            fprintf(fid, '  double *sv = data->sv;\n');
            fprintf(fid, '  double *dvdx = data->dvdx;\n');
            fprintf(fid, '  double *dvdu = data->dvdu;\n');
            fprintf(fid, '  double *dvdp = data->dvdp;\n');
            fprintf(fid, '  double *su = data->su;\n');
            % 	fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
            fprintf(fid, '  double *sx_tmp = N_VGetArrayPointer(sx);\n');
            fprintf(fid, '  double *sxdot_tmp = N_VGetArrayPointer(sxdot);\n');
            
            % Equations
            fprintf(fid, '  switch (ip) {\n');
            for j2=1:size(condition.sym.fsx,2)
                fprintf(fid, '  case %i: {\n', j2-1);
                if(timedebug)
                    fprintf(fid, '  printf("%%g \\t fsx%i\\n", t);\n', j2-1);
                end
                if(j2==1)
                    fprintf(fid, '  fsu_%s(data, t);\n', condition.fkt);
                    fprintf(fid, '  dvdx_%s(t, x, data);\n', condition.fkt);
                    fprintf(fid, '  dvdu_%s(t, x, data);\n', condition.fkt);
                    fprintf(fid, '  dvdp_%s(t, x, data);\n', condition.fkt);
                end
                writeCcode(fid, condition, 'fsv', j2);
                writeCcode(fid, condition, 'fsx', j2);
                fprintf(fid, '  } break;\n\n');
            end
            fprintf(fid, '  }\n');
            fprintf(fid, '  for (is=0; is<%i; is++) {\n', length(model.xs));
            fprintf(fid, '    if(mxIsNaN(sxdot_tmp[is])) sxdot_tmp[is] = 0.0;\n');
            fprintf(fid, '  }\n');
        end
    end
    fprintf(fid, '\n  return(0);\n}\n\n\n');
end


% write fsx0
fprintf(fid, ' void fsx0_%s(int ip, N_Vector sx0, void *user_data)\n{\n', condition.fkt);
if(~isempty(model.xs))
    if(config.useSensis)
        fprintf(fid, '  UserData data = (UserData) user_data;\n');
        fprintf(fid, '  double *p = data->p;\n');
        fprintf(fid, '  double *u = data->u;\n');
        fprintf(fid, '  double *sx0_tmp = N_VGetArrayPointer(sx0);\n');
        
        % Equations
        fprintf(fid, '  switch (ip) {\n');
        for j2=1:size(condition.sym.fsx0,2)
            fprintf(fid, '  case %i: {\n', j2-1);
            writeCcode(fid, condition, 'fsx0', j2);
            fprintf(fid, '  } break;\n\n');
        end
        fprintf(fid, '  }\n');
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');


% write dfxdp
fprintf(fid, ' void dfxdp_%s(realtype t, N_Vector x, double *dfxdp, void *user_data)\n{\n', condition.fkt);
if(timedebug) 
	fprintf(fid, '  printf("%%g \\t dfxdp\\n", t);\n');
end
if(~isempty(model.xs))
    if(config.useSensis)
        if(~isempty(condition.sym.dfxdp))
            fprintf(fid, '  int is;\n');
            fprintf(fid, '  UserData data = (UserData) user_data;\n');
            fprintf(fid, '  double *p = data->p;\n');
            fprintf(fid, '  double *u = data->u;\n');
            fprintf(fid, '  double *dvdp = data->dvdp;\n');
            fprintf(fid, '  double *dvdx = data->dvdx;\n');
            fprintf(fid, '  double *dvdu = data->dvdu;\n');
            fprintf(fid, '  double *x_tmp = N_VGetArrayPointer(x);\n');
            writeCcode(fid, condition, 'dfxdp');
            fprintf(fid, '  for (is=0; is<%i; is++) {\n', numel(condition.sym.dfxdp));
            fprintf(fid, '    if(mxIsNaN(dfxdp[is])) dfxdp[is] = 0.0;\n');
            fprintf(fid, '  }\n');
        end
    end
end
fprintf(fid, '\n  return;\n}\n\n\n');


% write OBS headers
function arWriteHFilesOBS(fid, data)

fprintf(fid, '#ifndef _MY_%s\n', data.fkt);
fprintf(fid, '#define _MY_%s\n\n', data.fkt);

fprintf(fid, '#include <cvodes/cvodes.h>\n'); 
fprintf(fid, '#include <cvodes/cvodes_dense.h>\n');
fprintf(fid, '#include <nvector/nvector_serial.h>\n');
fprintf(fid, '#include <sundials/sundials_types.h>\n'); 
fprintf(fid, '#include <sundials/sundials_math.h>\n');  
fprintf(fid, '#include <udata.h>\n');
fprintf(fid, '#include <math.h>\n');
fprintf(fid, '#include <mex.h>\n');
fprintf(fid, '#include <arInputFunctionsC.h>\n');
fprintf(fid,'\n\n\n');

fprintf(fid, ' void fy_%s(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int iruns, double *y, double *p, double *u, double *x);\n', data.fkt);
fprintf(fid, ' void fystd_%s(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x);\n', data.fkt);
fprintf(fid, ' void fsy_%s(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *su, double *sx);\n', data.fkt);
fprintf(fid, ' void fsystd_%s(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *sy, double *su, double *sx);\n\n', data.fkt);

fprintf(fid, '#endif /* _MY_%s */\n', data.fkt);
fprintf(fid,'\n\n\n');


% OBS
function arWriteCFilesOBS(fid, config, data)

fprintf(fid, '#include "%s.h"\n',  data.fkt);
fprintf(fid, '#include <cvodes/cvodes.h>\n');    
fprintf(fid, '#include <cvodes/cvodes_dense.h>\n');
fprintf(fid, '#include <nvector/nvector_serial.h>\n');
fprintf(fid, '#include <sundials/sundials_types.h>\n'); 
fprintf(fid, '#include <sundials/sundials_math.h>\n');  
fprintf(fid, '#include <udata.h>\n');
fprintf(fid, '#include <math.h>\n');
fprintf(fid, '#include <mex.h>\n');
fprintf(fid, '#include <arInputFunctionsC.h>\n');
fprintf(fid,'\n\n\n');

% write y
fprintf(fid, ' void fy_%s(double t, int nt, int it, int ntlink, int itlink, int ny, int nx, int iruns, double *y, double *p, double *u, double *x){\n', data.fkt);
writeCcode(fid, data, 'fy');
fprintf(fid, '\n  return;\n}\n\n\n');

% write ystd
fprintf(fid, ' void fystd_%s(double t, int nt, int it, int ntlink, int itlink, double *ystd, double *y, double *p, double *u, double *x){\n', data.fkt);
writeCcode(fid, data, 'fystd');
fprintf(fid, '\n  return;\n}\n\n\n');

% write sy
fprintf(fid, ' void fsy_%s(double t, int nt, int it, int ntlink, int itlink, double *sy, double *p, double *u, double *x, double *su, double *sx){\n', data.fkt);
if(config.useSensis)
    writeCcode(fid, data, 'fsy');
end
fprintf(fid, '\n  return;\n}\n\n\n');

% write systd
fprintf(fid, ' void fsystd_%s(double t, int nt, int it, int ntlink, int itlink, double *systd, double *p, double *y, double *u, double *x, double *sy, double *su, double *sx){\n', data.fkt);
if(config.useSensis)
    writeCcode(fid, data, 'fsystd');
end
fprintf(fid, '\n  return;\n}\n\n\n');


% write C code
function writeCcode(fid, cond_data, svar, ip)

if(strcmp(svar,'fv'))
    cstr = ccode(cond_data.sym.fv(:));
    cvar =  'data->v';
elseif(strcmp(svar,'dvdx'))
    cstr = ccode(cond_data.sym.dfvdx(:));
    cvar =  'data->dvdx';
elseif(strcmp(svar,'dvdu'))
    cstr = ccode(cond_data.sym.dfvdu(:));
    cvar =  'data->dvdu';
elseif(strcmp(svar,'dvdp'))
    cstr = ccode(cond_data.sym.dfvdp(:));
    cvar =  'data->dvdp';
elseif(strcmp(svar,'fx'))
    cstr = ccode(cond_data.sym.fx(:));
    for j=find(cond_data.sym.fx(:)' == 0)
        cstr = [cstr sprintf('\n  T[%i][0] = 0.0;',j-1)]; %#ok<AGROW>
    end
    cvar =  'xdot_tmp';
elseif(strcmp(svar,'fx0'))
    cstr = ccode(cond_data.sym.fpx0(:));
    cvar =  'x0_tmp';
elseif(strcmp(svar,'dfxdx'))
    cstr = ccode(cond_data.sym.dfxdx(:));
    for j=find(cond_data.sym.dfxdx(:)' == 0)
        cstr = [cstr sprintf('\n  T[%i][0] = 0.0;',j-1)]; %#ok<AGROW>
    end
    cvar =  'J->data';
elseif(strcmp(svar,'fsv'))
    cstr = ccode(cond_data.sym.fsv(:,ip));
    cvar =  'sv';
elseif(strcmp(svar,'fsx'))
    cstr = ccode(cond_data.sym.fsx(:,ip));
    for j=find(cond_data.sym.fsx(:,ip)' == 0)
        cstr = [cstr sprintf('\n  T[%i][0] = 0.0;',j-1)]; %#ok<AGROW>
    end
    cvar =  'sxdot_tmp';
elseif(strcmp(svar,'fsx0'))
    cstr = ccode(cond_data.sym.fsx0(:,ip));
    cvar =  'sx0_tmp';
elseif(strcmp(svar,'fu'))
    cstr = ccode(cond_data.sym.fu(:));
    cvar =  'data->u';
elseif(strcmp(svar,'fsu'))
    cstr = ccode(cond_data.sym.dfudp(:));
    cvar =  'data->su';
elseif(strcmp(svar,'fy'))
    cstr = ccode(cond_data.sym.fy(:));
    cvar =  'y';
elseif(strcmp(svar,'fystd'))
    cstr = ccode(cond_data.sym.fystd(:));
    cvar =  'ystd';
elseif(strcmp(svar,'fsy'))
    cstr = ccode(cond_data.sym.fsy(:));
    cvar =  'sy';
elseif(strcmp(svar,'fsystd'))
    cstr = ccode(cond_data.sym.fsystd(:));
    cvar =  'systd';
elseif(strcmp(svar,'dfxdp'))
    cstr = ccode(cond_data.sym.dfxdp(:));
    cvar =  'dfxdp';
end

cstr = strrep(cstr, 't0', [cvar '[0]']);
cstr = strrep(cstr, '][0]', ']');
cstr = strrep(cstr, 'T', cvar);

% % debug
% fprintf('\n\n');
% if(config.isMaple)
%     for j=1:length(cstr)
%         fprintf('%s\n', cstr{j});
%     end
% else
%     fprintf('%s', cstr);
% end
% fprintf('\n');

if(~(length(cstr)==1 && isempty(cstr{1})))
	if(strcmp(svar,'fy') || strcmp(svar,'fystd') || strcmp(svar,'fsy') || strcmp(svar,'fsystd'))
        if(strcmp(svar,'fy'))
            cstr = strrep(cstr, 'x[', 'x[nx*ntlink*iruns+itlink+ntlink*');
            cstr = strrep(cstr, 'y[', 'y[ny*nt*iruns+it+nt*');
        else
            cstr = strrep(cstr, 'x[', 'x[itlink+ntlink*');
            cstr = strrep(cstr, 'y[', 'y[it+nt*');
        end
		cstr = strrep(cstr, 'u[', 'u[itlink+ntlink*');
		cstr = strrep(cstr, 'ystd[', 'ystd[it+nt*');
	else
		cstr = strrep(cstr, 'x[', 'x_tmp[');
		cstr = strrep(cstr, 'dvdx_tmp', 'dvdx');
	end
end

fprintf(fid, '%s\n', cstr);

% % debug
% fprintf('\n\n%s\n', cstr);











