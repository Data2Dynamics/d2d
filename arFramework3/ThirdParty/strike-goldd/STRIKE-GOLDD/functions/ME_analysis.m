%=========================================================================%
%========ME_analysis: Convert a model into Multi-Experiment form==========%

function ME_analysis(modelname,opts)

%============================Load model===================================%
load(modelname); %#ok<*LOAD>

%======================Dimensions of the problem==========================%
% Number of states:
n=numel(x); %#ok<*NODEF>
% Number of outputs:
m=numel(h);
% Number of known inputs:
if exist('u','var')
    nu=numel(u);
else
    u=[];
    nu=0;
end
% Number of unknown inputs:
if exist('w','var')
    nw=numel(w);
else
    w=[];
    nw=0;
end
%========================Initialize variables=============================%

me_x   = sym(zeros(n*opts.multiexp_numexp,1));
me_h   = sym(zeros(m*opts.multiexp_numexp,1));
me_f   = sym(zeros(n*opts.multiexp_numexp,1));
me_u   = sym(zeros(nu*opts.multiexp_numexp,1));
me_w   = sym(zeros(nw*opts.multiexp_numexp,1));

%========================Multi-experiment model===========================%

variables=[x;w;u];

for i=1:opts.multiexp_numexp
    
    num_exp=strcat('Exp',num2str(i));
    
    for ind_u=1:nu,me_u(ind_u+nu*(i-1)) = sym(strcat(char(u(ind_u)),sprintf(num_exp)));end
    for ind_w=1:nw,me_w(ind_w+nw*(i-1)) = sym(strcat(char(w(ind_w)),sprintf(num_exp)));end
    for ind_x=1:n,me_x(ind_x+n*(i-1))= sym(strcat(char(x(ind_x)),sprintf(num_exp)));end

    aug_variables=[me_x(1+n*(i-1):n*i);me_w(1+nw*(i-1):nw*i);me_u(1+nu*(i-1):nu*i)];
    
    for ind_f=1:n,me_f(ind_f+n*(i-1))= subs(f(ind_f),variables,aug_variables);end
    for ind_h=1:m,me_h(ind_h+m*(i-1))= subs(h(ind_h),variables,aug_variables);end
end

% New variables: 
u=me_u;
w=me_w;
x=me_x;
% Multi-experiment dynamics:
f=me_f;
% Multi-experiment output:
h=me_h;

% Initial conditions:
if exist('ics','var')
    if opts.multiexp_user_ics == 0
        % Replicate initial conditions:
        ics=repmat(ics,1,opts.multiexp_numexp);
        known_ics=repmat(known_ics,1,opts.multiexp_numexp);
    else
        ics=reshape(opts.multiexp_ics,1,[]);
        known_ics=reshape(opts.multiexp_known_ics,1,[]);
    end
else
    ics=[];
    known_ics=[];
end

if exist('p','var')==0
    p=[];
end

save(strcat(pwd,filesep,'models',filesep,strcat(modelname,'_',num2str(opts.multiexp_numexp)),'Exp'),'x','p','u','w','f','h','ics','known_ics')  
