
function [transf,nuevas_variables,allVar,z_v]=Lie_Symmetry(varargin)

    if nargin > 0
        modelname = varargin;
        load(modelname{1});
        [~,paths,opts,submodels,prev_ident_pars] = options();
        ansatz=opts.ansatz;
        pMax=opts.degree;
        tmax=opts.tmax;
        ode_n=opts.ode_n;
    else
        [modelname,paths,opts,submodels,prev_ident_pars] = options();
        load(modelname);
        ansatz=opts.ansatz;
        pMax=opts.degree;
        tmax=opts.tmax;
        ode_n=opts.ode_n;
    end

%clc
%%  SYMMETRY SEARCH ALGORITHM WITH INFINITESIMALS
%   Code for symmetry search algorithm with infinitesimals generators.
%   This is the procedure that is followed:
%   1. Reading data.
%   2. Creation of the infinitesimals polynomials.
%   3. Obtaining derivatives of the infinitesimal polynomials.
%   5. Obtaining derivatives of the numerator, denominator and ICS.
%   6. Construction of state polynomials, observable functions and ICS.
%   7. System resolution.
%   8. Calculate generators and transformations (if exist).
% 

%%   Reading data
if exist('u') && isempty(u)==0
    assume(u,'real');
end
%   Check if are columns or rows
if isrow(p)==1
    p=transpose(p);
end
if isrow(x)==1
   x=transpose(x); 
end
if exist('w','var') && isempty(w)==0
    if isrow(w)==1
        w=transpose(w);
    end
    assume(w,'real');
    allVar=[x;w;p];
    nw=length(w);
else
    allVar=[x;p];
    nw=0;
end
if exist('ics','var') && isempty(ics)==0
    nics=length(ics);
    ic=[];
    if isrow(ics)==1
         ics=transpose(ics); 
    end
    %   Numerator and denominator
    [num_ics,den_ics]=numden(sym(ics));
    ind=zeros(nics,1);
    %   Extract the coefficients
    for i=1:nics
        if (isempty(diff(ics(i)))==0 && diff(ics(i))==1)&& known_ics(i)==1
            assume(ics(i),'real');
            [cn,tn]=coeffs(num_ics(i),allVar);
            [cd,td]=coeffs(den_ics(i),allVar);
            [Lin,Lon] = ismember(cn,allVar);
            [Lid,Lod] = ismember(cd,allVar);
            l_Lin=length(Lin);
            l_Lid=length(Lid);
            for k=1:l_Lin
               if Lin(k)==0 && cn(k)~=1
                  allVar=[allVar;cn(k)]; 
                  ind(i)=i;
                  ic=[ic;cn(k)];
               end
            end
            for k=1:l_Lid
                if Lid(k)==0 && cd(k)~=1
                   allVar=[allVar;cd(k)]; 
                end
            end
        end
    end
    %   Obtaining the variables within ICS
    ind=nonzeros(ind);
    nics=length(ind);
    ics_n=subs(ics,x(ind),ic);
    %   Replace recursively
    cont=0;
    while cont<3
        ics_n=subs(ics_n,x,ics_n);
        cont=cont+1;
    end
    [num_ics,den_ics]=numden(ics_n);
else
    nics=0;
end


%   Creation variables and constants
l_x=length(x);              %   States vector length
l_p=length(p)+nics;         %   Parameters vector length
l_f=length(f);              %   f function vector length
l_h=length(h);              %   h function vector length
l_allVar=length(allVar);    %   allVar length
assume(allVar,'real');

%	Start of computing time
tic

%%  ANSATZ
%   Creation of infinitesimal polynomials
if nics~=0
    m_aux=[1;p;ic];
else
    m_aux=[1;p];
end
if (ansatz==1)
    [infi,rs_k] = create_uni(allVar,l_allVar,pMax);
elseif (ansatz==2)
    [infi,rs_k] = create_par(allVar,m_aux,l_p,l_x,p,pMax,nw);
else
    [infi,rs_k] = create_multi(allVar,m_aux,l_p,l_x,p,pMax,nw,x);
end
l_rs=length(rs_k);

%%  ANSATZ DERIVATIVES
diffini=[];
if (ansatz==3)
    for i=1:length(infi)
        diffi=[];
        for j=1:l_allVar
            diffi=[diffi,diff(infi(i),allVar(j))];
        end
        diffini=[diffini;diffi];
    end
else
    for i=1:l_allVar
        diffini=[diffini,diff(infi(i),allVar(i))];
    end
end
%% NUMERATOR AND DENOMINATOR
[num_f,den_f]= numden(f);
[num_h,den_h]= numden(h);

%% DERIVATIVES OF THE NUMERATOR, DENOMINATOR AND ICS
diff_num_f = [];
diff_num_h = [];
diff_den_f = [];
diff_den_h = [];
if isempty(ics)==0
    diff_num_ics=[];
    diff_den_ics=[];
    for i=1:l_allVar
       diff_num_f = [diff_num_f,diff(num_f,allVar(i))];
       diff_den_f = [diff_den_f,diff(den_f,allVar(i))];
       diff_num_h = [diff_num_h,diff(num_h,allVar(i))];
       diff_den_h = [diff_den_h,diff(den_h,allVar(i))];
       diff_num_ics=[diff_num_ics,diff(num_ics,allVar(i))];
       diff_den_ics=[diff_den_ics,diff(den_ics,allVar(i))];
    end
else
    for i=1:l_allVar
       diff_num_f = [diff_num_f,diff(num_f,allVar(i))];
       diff_den_f = [diff_den_f,diff(den_f,allVar(i))];
       diff_num_h = [diff_num_h,diff(num_h,allVar(i))];
       diff_den_h = [diff_den_h,diff(den_h,allVar(i))];
    end
end    

%% CONSTRUCTION OF POLYNOMIES AND SYSTEM
l_in=length(infi);
%   POL. STATES
if (ansatz==3)
    [A_sta] = pol_sta_2(allVar,den_f,diffini,diff_den_f,diff_num_f,...
                        infi,l_f,l_in,num_f,rs_k);
else
    [A_sta] = pol_sta_1(allVar,den_f,diffini,diff_den_f,diff_num_f,...
                        infi,l_f,l_in,num_f,rs_k);
end

%   POL. OBS.
[A_obs] = pol_obs(allVar,den_h,diff_den_h,diff_num_h,...
                        infi,l_h,l_in,num_h,rs_k);                  
                    
%   ICS
if isempty(ics)==0 
    [A_ics] = pol_ics(allVar,den_ics,diff_den_ics,diff_num_ics,ics_n,...
                infi,known_ics,l_in,l_x,num_ics,rs_k,x);
    %   SYSTEM
    A=[A_sta;A_obs;A_ics];
else
    %   SYSTEM
    A=[A_sta;A_obs];
end

%%  SYSTEM RESOLUTION
V=null(A);
si2=size(V);
ind=si2(2);

%%  TRANSFORMATIONS (IF THEY EXIST)
if ind~=0
    % There is/are transformation(s) with this Ansatz and degree
    [transf,nuevas_variables,z_v] = new_var(V,ind,infi,rs_k,allVar,l_allVar,tmax,ode_n);    
else
    fprintf('\n--------------------------------------------------\n');
    fprintf('>>> There are no transformations with these degree and/or Ansatz \n')
    fprintf('\n--------------------------------------------------\n');
    transf=[];
    nuevas_variables=[];
    allVar=[];
    z_v=[];
    return
end
toc

end