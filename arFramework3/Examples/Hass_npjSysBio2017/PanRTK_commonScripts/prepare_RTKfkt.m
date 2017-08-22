%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   prepare_RTKfkt                          %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Part of estimation of initial concentrations for signaling model
%  prepare equation system with desired receptor/surface level
%
%
function prepare_RTKfkt(RNAseq)
    global ar;
    if(length(RNAseq)~=length(ar.intials))
        error('provide RNAseq vector of length of RTKs');
    end
    ar.initials_tmp = ar.intials;
    for i=1:length(ar.intials)
        ar.initials_tmp{i} = sprintf('%s - %d',ar.intials{i},10^RNAseq(i));
    end
end