function [new_var] =sol(allVar, transf) 
%% FUNCTION TO COMPUTE THE ODE SOLUTION OF THE INFINITESIMAL GENERATOR
    k=find(transf);
    ttransf=transf(k);
    l=length(k);
    syms aa(epsilon) [1 l]
    aV=allVar(k)';
    % Create equations
    eqqs=[];
    var_a=[];
    t_ch=subs(ttransf, aV,aa);
    eqs=diff(aa(epsilon),epsilon)==t_ch;
    cond=aa(0)==aV;
    eqs=transpose(eqs);
    cond=transpose(cond);
    S=dsolve(eqs,cond);
    snv=[];
    AA=sym('aa', [l 1]);
    AA=string(AA);
    S=orderfields(S,AA);
    v=fieldnames(S);
    for j=1:length(v)
        myVar=S.(v{j});
        snv=[snv,myVar];
    end
    new_var=allVar;
    new_var(k)=snv;
end

