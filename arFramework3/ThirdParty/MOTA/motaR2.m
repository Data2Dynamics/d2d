function r2=motaR2(R,Ps)
% r2=motaR2(R,Ps)
% motaR2 expects R and Ps to be ordered along the comlumns, i.e., different
% predictor variables a ordered along the rows.
% R  = Response
% Ps = Predictors

%% mean response
    meanR=mean(R);
%% SS-Residual
   if length(Ps(1,:))~=1
      SS_res=sum((R-sum(Ps,2)).^2);
   else
   
      SS_res=sum((R-Ps).^2);
   end
    
%% SS-Total 
   SS_tot=sum((R-meanR).^2);
%% r2
   r2=1-(SS_res/SS_tot);
end