%--------------------------------------------------------------------------
% Function that builds the generalized observability-identifiability matrix
% of the model specified in the 'options' script. It takes as input the
% number of Lie derivatives, provided by the user. The resulting array is
% stored in a MAT file.
%--------------------------------------------------------------------------

function build_OI_ext(nd)

%==========================================================================
% Read options, add folders to path:
global x f p u w wlvector
[modelname,paths,opts,submodels,prev_ident_pars] = options; %#ok<*ASGLU>
addpath(genpath(paths.models));
addpath(genpath(paths.results));
addpath(genpath(paths.functions));

%==========================================================================
% Load model:
load(modelname)
fprintf('\n Analyzing identifiability of %s ... \n', modelname);
    
% ==========================================================================
% Remove parameters that have already been classified as identifiable:
if exist('prev_ident_pars','var') == 1
    for np=1:numel(prev_ident_pars)
        [~, original_index] = ismember(prev_ident_pars(np),p); 
        p(original_index)=[];                
    end
end

%==========================================================================
% Dimensions of the problem: 
m    = numel(h);                  % number of outputs
n    = numel(x);                  % number of states
q    = numel(p);                  % number of unknown parameters
if exist('w','var')
    nw = numel(w); % number of unknown inputs
else 
    nw = 0;
    w = [];
end
if exist('u','var')
    nu = numel(u); % number of unknown inputs
else 
    nu = 0;
    u = [];
end
fprintf('\n >>> The model contains:\n %d states:\n %s',n,char(x));
fprintf('\n %d outputs:\n %s',m,char(h));
fprintf('\n %d known inputs:\n %s',nu,char(u));
fprintf('\n %d unknown inputs:\n %s',nw,char(w));
fprintf('\n %d parameters:\n %s',q,char(p));

%==========================================================================
fprintf('\n >>> Building Identifiability-Observability matrix with %d Lie derivatives...',nd);
tic

%======================================================================
% Input derivatives:

%- Create array of known inputs and set certain derivatives to zero:
if numel(u)>0
    for ind_u=1:numel(u) % create array of derivatives of the inputs
        input_der(ind_u,:) = [u(ind_u),sym(strcat(char(u(ind_u)),sprintf('_d')),[1 nd-1])];  %[1 max(opts.nnzDerU)]
        input_der(ind_u,(opts.nnzDerU(ind_u)+2):end)=0;
    end          
else
    input_der = [];
end
syms zero_input_der_dummy_name

%- Create array of unknown inputs and set certain derivatives to zero:
if numel(w)>0
    for ind_w=1:numel(w) 
        w_der(ind_w,:) = [w(ind_w),sym(strcat(char(w(ind_w)),sprintf('_d')),[1 nd+1])];  % [1 max(opts.nnzDerW)]
        w_der(ind_w,(opts.nnzDerW(ind_w)+2):end)=0;
    end 
    wlvector     = reshape(w_der(:,1:end-1),[],1);  % reshape array as a column vector
    wlvector_dot = reshape(w_der(:,2:end),[],1);    % vector of derivatives of the unknown inputs 
    %-- Include as states only nonzero inputs/derivatives:
    [nzi,nzj,nz_wlvec] = find(wlvector); %#ok<ASGLU>
    wlvector     = nz_wlvec;   
    wlvector_dot = wlvector_dot(nzi); 
else
    wlvector      = [];
    wlvector_dot  = [];
end  

%======================================================================
% Augment state vector, dynamics:
xaug = [x;p;wlvector];                      
faug = [f;zeros(numel(p),1);wlvector_dot]; 

%======================================================================
% Build Oi:
onx        = zeros(m*(1+nd),n+q+numel(wlvector)); 
onx        = sym(onx);
onx(1:m,:) = jacobian(h,xaug); % 1st block
totaltime  = toc;
ind        = 0; % Lie derivative index (k)
lasttime   = 0; 
%----------------------------------------------------------------------
past_Lie   = h;        
extra_term = 0;       
while ind < nd && lasttime < opts.maxLietime % 2nd and subsequent blocks
    tic
    Lieh = onx((ind*m+1):(ind+1)*m,:)*faug;
    if ind>0 
        if numel(u) > 0
            for i=1:ind 
                if i < size(input_der,2) 
                    lo_u_der   = input_der(:,i);   
                    hi_u_der   = input_der(:,i+1);
                    lo_u_der   = subs(lo_u_der,0,zero_input_der_dummy_name);
                    extra_term = extra_term + jacobian(past_Lie,lo_u_der)*hi_u_der;
                end                        
            end
        end
    end       
    ext_Lie = Lieh + extra_term;
    past_Lie = ext_Lie;         
    onx(((ind+1)*m+1):((ind+2)*m),:) = jacobian(ext_Lie,xaug); 
    lasttime = toc;
    totaltime = totaltime + lasttime;
    ind = ind+1;
    obsidentmatrix = sprintf('obs_ident_matrix_%s_%d_Lie_deriv',modelname,ind);
    obsidentfile = strcat(pwd,filesep,'results',filesep,obsidentmatrix);
    save(obsidentfile);
    fprintf('\n >>> Lie derivative %d calculated in %d seconds',ind,totaltime);
end
fprintf('\n\n');
end    