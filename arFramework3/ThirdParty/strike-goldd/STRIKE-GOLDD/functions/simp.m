function [sim_vec] = simp(allVar,vec,eps_var)
%% FUNCTION SIMPLIFY
%
%   Function to find some expressions that may have not been simplified
%   during the computation of new variables
%
%   INPUT:      allVar --> All variables
%               vec --> vector with final tranformations
%               eps_var --> value of exponential(epsilon)
%
%   OUTPUT:
%               sim_vec --> simplified final transformations vector
%%
    syms epsilon
    ch=children(vec);
    [~,ncol]=size(ch);
    if iscell(ch)==0 % For R2020b
        ch=num2cell(ch);
    end
    for i=1:ncol
       [~,ncol_2]=size(ch{1,i});
       if iscell(ch{1,i})==1 %For R2020b
           aux_v=[ch{1,i}{:}];
       else
           aux_v=[ch{1,i}(:)];
       end
       for j=1:ncol_2
          [n,d]=numden(aux_v(j));
          if n==1 %% Only denominator 
              [c,p]=coeffs(d,allVar);
              if (numel(unique(abs(p)))==1 && abs(c)~=1) % Some expression contains exp(ep)
                  eq=simplify(log(c),'IgnoreAnalyticConstraints',true);
                  [cc,~]=coeffs(eq,epsilon);
                  vec(i)=subs(vec(i),exp(cc*epsilon),eps_var^(cc));
              end
          elseif d==1
              [c,p]=coeffs(n,allVar);
              if (numel(unique(abs(p)))==1 && abs(c)~=1) % Some expression contains exp(ep)
                  eq=simplify(log(c),'IgnoreAnalyticConstraints',true);
                  [cc,~]=coeffs(eq,epsilon);
                  vec(i)=subs(vec(i),exp(cc*epsilon),eps_var^(cc));
              end
          else
              [c_n,p_n]=coeffs(n,allVar);
              [c_d,p_d]=coeffs(d,allVar);
              if numel(unique(abs(p_n)))==1 && abs(c_n)~=1 % Some expression contains exp(ep)
                  eq=simplify(log(c_n),'IgnoreAnalyticConstraints',true);
                  [cc_n,~]=coeffs(eq,epsilon);
                  vec(i)=subs(vec(i),exp(cc_n*epsilon),eps_var^(cc_n));
              elseif numel(unique(abs(p_d)))==1 && abs(c_d)~=1 
                  eq=simplify(log(c_d),'IgnoreAnalyticConstraints',true);
                  [cc_d,~]=coeffs(eq,epsilon);
                  vec(i)=subs(vec(i),exp(cc_d*epsilon),eps_var^(cc_d));
              end
          end
       end
    end
    sim_vec=vec;
end

