%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   crawl_viability                         %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Go through all cell lines in CCLE with their respective receptor
%  expression values, calculates their model response and fills matrices
%  that fit to those obtained in cell viability measurements

function crawl_viability(CCLE, norm,do_qFACS)
global ar;

if(~exist('norm','var') || isempty(norm))
   norm = 1;
end
if(~exist('do_qFACS','var') || isempty(do_qFACS))
   if(size(CCLE,1)==5)
       do_qFACS = 0;
   else
       do_qFACS=1;
   end
end
if(size(CCLE,1)==5)
    cell_line_expr = CCLE(1:5,:);
else
    cell_line_expr = CCLE([1:4 6],:);
end

init_vec={'EGFR','ErbB2','ErbB3','IGF1R','Met'};
RTK_level=[];
for i=1:length(init_vec)
    RTK_level=[RTK_level strmatch(strcat('init_',init_vec{i}),ar.pLabel,'exact')];
end

for i=1:size(cell_line_expr,2)
    %Get init parameters for the corresponding RNAseq
    inits=get_initals(cell_line_expr(:,i),do_qFACS);
    
    ar.p(RTK_level)=inits;
    if(sum(inits==-1)==length(inits))
        continue;
    end
    %get the model response for the MM drugs
    get_ABDiffs()
    get_model_response('ALL',i,norm);
     
end

end