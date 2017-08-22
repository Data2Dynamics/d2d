%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   calc_RTK                                %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Part of estimation of initial concentrations for signaling model
%  Solve equation system to obtain intial receptor concentrations that have
%  to be set in model to obtain desired receptor/surface levels
%

function out=calc_RTK(x,init_vec)
global ar;
if(~exist('x','var') || isempty(x))
    error('I dont get it')
end
if(length(x) ~= length(init_vec))
    error('Wrong number of RNAseq specified');
end
out = NaN(size(ar.intials));
for i=1:length(ar.intials)
    fkt_tmp = ar.initials_tmp{i};
    %symvar(fkt_tmp);
    for j=1:length(x)
        fkt_tmp=strrep(fkt_tmp,init_vec{j},num2str(10^x(j)));
        
    end
    out(i)=str2num(fkt_tmp);
end
