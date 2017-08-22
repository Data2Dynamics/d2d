%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   get_initials                            %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Main script for estimation of initial concentrations for signaling model
%  Take receptor expression or surface levels as input and solve
%  observations of receptors on surface to obtain init_XX for model
%

function inits=get_initals(RNAseq, do_qFACS)
global ar

options = optimoptions('fsolve','Algorithm','levenberg-marquardt','DiffMinChange',1.e-4,'Display','off');
if(~isfield(ar,'intials') || sum(isempty(ar.intials))>0)
    get_RTKlevel;
end

if(~exist('do_qFACS','var') || isempty(do_qFACS))
   do_qFACS = 1; 
end

init_vec={'EGFR','ErbB2','ErbB3','IGF1R','Met'};
is_there = ismember(strcat('init_',init_vec),ar.pLabel);
init_vec = init_vec(is_there);

%Calc receptor surface from RNAseq
if(do_qFACS)
    qFACS(1) = log10(sinh(asinh(RNAseq(1))*0.964 + 8.9229)/10000);
    qFACS(2) = log10(sinh(asinh(RNAseq(2))*0.8282 + 8.0359)/10000);
    qFACS(3) = log10(sinh(asinh(RNAseq(3))*0.3445 + 9.1685)/10000);
    if(is_there(4)==1)
        qFACS(4) = log10(sinh(asinh(RNAseq(4))*0.4757 + 9.47)/10000);
    elseif(is_there(5)==1)
        qFACS(4) = log10(sinh(asinh(RNAseq(4))*0.6418 + 8.3577)/10000);
    end
    if(is_there(4)==1 && is_there(5)==1)
        qFACS(5) = log10(sinh(asinh(RNAseq(5))*0.6418 + 8.3577)/10000);
    end
else
   qFACS = RNAseq; 
end
prepare_RTKfkt(qFACS);
try
    inits = fsolve(@calc_RTK,ones(size(qFACS))*-1,options,init_vec);
catch
    sprintf('something went wrong in getting inits, retrying')
    try
        inits = fsolve(@calc_RTK,ones(size(qFACS))*1,options,init_vec);
    catch
        inits = ones(size(qFACS))*-1;
    end
end