% This function compiles the NL2SOL code with the mex wrapper.
% It should be called to produce a mex file which is callable. Note that this requires
% both a C compiler as well as a FORTRAN compiler to be available on the system.
%
% NL2SOL is part of the PORT library. More specifically: 
%		Algorithm 573:  NL2SOL—An Adaptive Nonlinear Least-Squares 
%		Authors:        John Dennis, David Gay, Roy Welsch
%
% After compilation the function can be called (mostly) analogously to lsqnonlin.
%   e.g. [X,RESNORM,RESIDUAL,EXITFLAG,ITERATIONS,FEVALS,JACOBIAN] = 
%            mexnl2sol(@fun, x0, (lb), (ub), (opts), (printlevel))
%

function compileNL2SOL()

    % Compile NL2SOL
    cpath   = mfilename('fullpath');
    loc     = strfind( fliplr(cpath), '/');
    cpath   = cpath(1:end-loc+1);
    mask    = sprintf('%sPORT/%%s', cpath);
    
    F = {   'dn2gb.f', 'dn2g.f', 'dn2f.f', 'dn2fb.f', 'dh2rfg.f'...
            'Mach/d1mach.f', 'seterr.f', 'drn2gb.f', 'stopx.f', 'dd7mlp.f', ...
            'drn2g.f', 'dv7scp.f', 'dn2rdp.f', 'dv7dfl.f', 'dh2rfa.f', ...
            'divset.f', 'dr7mdc.f', 's88fmt.f', 'i7shft.f', 'dl7msb.f', ...
            'dq7rsh.f', 'dv7vmp.f', 'dl7sqr.f', 'ds7ipr.f', 'dg7qsb.f'...
            'i1mach.f', 'e9rint.f', 'i8save.f', 'eprint.f', 'ds7bqn.f', ...
            'sdump.f', 'ds7dmp.f', 'drldst.f', 'da7sst.f', 'ds7lvm.f', ...
            'ds7lup.f', 'df7dhb.f', 'df7hes.f', 'dl7srt.f', 'dg7itb.f', ...
            'dv7cpy.f', 'dq7apl.f', 'dl7tvm.f', 'dv2nrm.f', 'dr7tvm.f', ...
            'dd7tpr.f', 'dd7upd.f', 'dl7svn.f', 'dl7mst.f', 'dl7vml.f', ...
            'dl7svx.f', 'dg7lit.f', 'dv7scl.f', 'dq7rad.f', 'dparck.f', ...
            'dl7itv.f', 'dv7ipr.f', 'do7prd.f', 'dn2lrd.f', 'i7pnvr.f', ...
            'i7mdcn.f', 'dl7ivm.f', 'dv2axy.f', 'dn2cvp.f', 'fdump.f', ...
            'i0tk00.f', 'iceil.f', 'iflr.f', 'i10wid.f', 'frmatr.f', ...
            'frmatd.f', 'frmati.f', 'u9dmp.f', 'a9rntl.f', 'a9rnti.f', ...
            'a9rntr.f', 'a9rntd.f', 'a9rntc.f', 'stkdmp.f', 'dg7qts.f', ...
            'dl7nvr.f', 'dl7tsq.f', 'i7copy.f', 'dc7vfn.f', 'ditsum.f', 'dv7shf.f' ...
            };

    fprintf( '\nCompiling NL2SOL . . . ' );
    files = cellfun(@(fn)sprintf(mask,fn), F(:), 'UniformOutput', false);
    outFiles = cellfun(@(fn)strrep(getFileName(fn), '.f', '.o'), F, 'UniformOutput', false);
    mex( '-c', '-largeArrayDims', '-lmwblas', '-lmwlapack', files{:} );
    
    try
        cc = mex.getCompilerConfigurations('fortran');
    catch
        error( 'No Fortran Compiler installed' );
    end
    
    if ( strcmp( cc.Name, 'gfortran' ) )
        mex( '-largeArrayDims', '-lgfortran', '-lmwblas', '-lmwlapack', sprintf('%s/mexnl2sol.c', cpath), outFiles{:}, '-outdir', cpath);
    else
        mex( '-largeArrayDims', '-lmwblas', '-lmwlapack', sprintf('%s/mexnl2sol.c', cpath), outFiles{:}, '-outdir', cpath);
    end
    
    delete(outFiles{:});
    fprintf( '[ OK ]\n');
end

function filename = getFileName( fullname )
    [~, filename, ext] = fileparts( fullname );
    filename = [filename ext];
end