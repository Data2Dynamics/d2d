%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   get_RTKlevel                            %%%
%%%   Helge Hass, 2016                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Part of estimation of initial concentrations for signaling model
%  Extract functions from model that calculate receptors/surface and
%  substitute with respective steady-state formulas, leaves initial
%  receptor concentrations as only degrees of freedom
%

function get_RTKlevel(get_total)
global ar;

if(~exist('get_total','var') || isempty(get_total))
    get_total = 0;
end

%Get pars for total RTK concentration
RTK_level=[];
init_vec={'EGFR','ErbB2','ErbB3','IGF1R','Met'};
for i=1:length(init_vec)
    RTK_level=[RTK_level strmatch(strcat('init_',init_vec{i}),ar.pLabel,'exact')];
end
ar.intials = cell(length(RTK_level),1);

phospho_vec={'init_EGFRi','init_ErbB2i','init_ErbB3i','init_IGF1Ri','init_pEGFRd', 'init_pEGFRi', 'init_pEGFRi_ph','init_pErbB12','init_pErbB12i','init_pErbB12i_ph','init_pErbB13','init_pErbB13i','init_pErbB13i_ph','init_pErbB2',...
    'init_pErbB2i','init_pErbB2i_ph','init_pErbB32','init_pErbB32i','init_pErbB32i_ph','init_pErbB3d','init_pErbB3i','init_pErbB3i_ph','init_pErbBi2','init_pErbBi2i','init_pErbBi2i_ph','init_pIGF1Rd','init_pIGF1Ri','init_pIGF1Ri_ph',...
    'init_Meti','init_pMetErbB3','init_pMetErbB3i','init_pMetErbB3i_ph','init_pMetd','init_pMetEGFR','init_pMeti_ph','init_pMeti','init_pMetEGFRi','init_pMetEGFRi_ph'};
p_vec=[];
for i=1:length(phospho_vec)
    p_vec = [p_vec strmatch(phospho_vec{i},ar.model(1).p,'exact')];
end

init_RTK = ar.model(1).p(p_vec);
RTK_fkt = ar.model(1).fp(p_vec);

%Get right fz
if(~get_total)
    fzs = find(~cellfun(@isempty,strfind(ar.model(1).z,'facs_')));
else
    fzs = [find(~cellfun(@isempty,strfind(ar.model(1).z,'total_'))) find(~cellfun(@isempty,strfind(ar.model(1).z,'_igf1r'))) find(~cellfun(@isempty,strfind(ar.model(1).z,'_met')))];
end

for i=1:length(fzs)
    expr_tmp = ar.model(1).fz{fzs(i)};
    
    for j=1:length(init_RTK)
        for k=1:length(ar.pLabel)
            if(sum(RTK_level==k)>0)
                continue
            end
            regPar = sprintf('(?<=\\W)%s(?=\\W)|^%s$|(?<=\\W)%s$|^%s(?=\\W)',ar.pLabel{k},ar.pLabel{k},ar.pLabel{k},ar.pLabel{k});
            if(ar.qLog10(k))
                RTK_fkt{j} = regexprep(RTK_fkt{j},regPar, num2str(10^ar.p(k)));
            else
                RTK_fkt{j} = regexprep(RTK_fkt{j},regPar, num2str(ar.p(k))); 
            end
        end
        get_Cond = ar.model(1).pcond(strmatch('is',ar.model(1).pcond));
        for k=1:length(get_Cond)
             regPar = sprintf('(?<=\\W)%s(?=\\W)|^%s$|(?<=\\W)%s$|^%s(?=\\W)',get_Cond{k},get_Cond{k},get_Cond{k},get_Cond{k});
            RTK_fkt{j} = regexprep(RTK_fkt{j},regPar, num2str(0));
        end      
        
        regPar = sprintf('(?<=\\W)%s(?=\\W)|^%s$|(?<=\\W)%s$|^%s(?=\\W)',strrep(init_RTK{j},'init_',''),strrep(init_RTK{j},'init_',''),strrep(init_RTK{j},'init_',''),strrep(init_RTK{j},'init_',''));
        expr_tmp = regexprep(expr_tmp,regPar, RTK_fkt{j});
    end 
    
    %part for only H322M replace
    reltos = strmatch('relto_',ar.model(1).data(1).pcond);
    for j=reltos'
        expr_tmp = strrep(expr_tmp,ar.model(1).data(1).pcond{j}, num2str(0));
    end
    
    expr_tmp = strrep(expr_tmp,'EGFR_EGF',num2str(0)); expr_tmp = strrep(expr_tmp,'EGFR_BTC',num2str(0)); expr_tmp = strrep(expr_tmp,'ErbB3_HRG',num2str(0)); expr_tmp = strrep(expr_tmp,'IGF1R_IGF1',num2str(0)); expr_tmp = strrep(expr_tmp,'Met_HGF',num2str(0));
    
    expr_tmp = strrep(expr_tmp,'init_','');
    
    ar.intials{i}=expr_tmp;
end