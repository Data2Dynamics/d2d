function r = create_heatmap_dr_lilly_epcam_var(ar2,interval,single,egf_in,min_val_met,max_val_met,min_val_ep,max_val_ep)

global ar
global result

ar = ar2;

if(~exist('min_val_met','var'))
    min_val_met = 100;
end
if(~exist('max_val_met','var'))
    max_val_met = 1000000;
end
if(~exist('min_val_ep','var'))
    min_val_ep = 100;
end
if(~exist('max_val_ep','var'))
    max_val_ep = 1000000;
end
if(~exist('interval','var'))
    interval = 60;
end
if(~exist('single','var'))
    single = 0;
end
if(~exist('egf_in','var'))
    egf_in = 0;
end


%x_doses = [300 100 33.333 11.111 3.7037 1.2346 0.41152 0.13717 0.045725];
x_doses = [1200 600 300 100 33.333 11.111 3.7037 1.2346 0.41152 0.13717 0.045725 0.01 0.001];
%x_doses = [300 100 33.333 11.111 3.7037 1.2346 0.41152 0.13717];
x_doses = log10(x_doses);
y_tmp(1:length(x_doses)) = 0;
ic90_req = 4.5;
ic50_rq = 22.5;

%no log!
arSetPars('init_egf',egf_in,0,0,-1,egf_in+1);   
%

arSetPars('init_egfr',log10((315000)),0,1,-5,7);
arSetPars('init_met',log10((71969.5)),0,1,-5,8);
arSetPars('init_epcam',log10((24045.75)),0,1,-5,8);

ar.model.qPlotYs(:) = 0;
ar.model.qPlotVs(:) = 0;

%read-out for ic90_doses  
ar.model.qPlotXs(:) = 0;
ar.model.qPlotXs(2) = 1;

arPlot
h = gcf;
axesObjs = get(h,'Children');
dataObjs = get(axesObjs,'Children');

for i = 1:length(ar.model.xNames)
    if(strcmp(ar.model.xNames(i),'pAkt'))
        pos_akt = length(axesObjs)-i;
    end
end

for k = 1:5
    ic_90_ro_x = dataObjs{pos_akt}(k).XData;
    ic_90_ro_y = dataObjs{pos_akt}(k).YData;
    if(length(ic_90_ro_x) == 830)
                ic90_req_10_min_val(k) = ic_90_ro_y(409);
    elseif(length(ic_90_ro_x) == 590)
                ic90_req_10_min_val(k) = ic_90_ro_y(192);
    else
                sprintf('New Simu_length( %i ), abording!', length(ic_90_ro_x))
                return;
    end
end

ic90_req_10_min_val = sort(ic90_req_10_min_val);
ic90_req_10_min_val = ic90_req_10_min_val(k) * 0.1;


sig_func_two = @(A,x)(A(4)+(A(1)-A(4))./((1+(x/A(3)).^A(2))));

[ic90_dose, ic50_dose] = deal(zeros(interval,interval));
die_ps_jc_fit_tmpy = zeros(interval,interval,length(ar.p));

[p_Akt_auc_mm131_jc_mm131_fit_tmpy, ys_ic_curve_mm131_jc_mm131_fit_tmpy, p_Akt_10_min_val_mm131_jc_mm131_fit_tmpy,p_Akt_30_min_val_mm131_jc_mm131_fit_tmpy,p_Akt_60_min_val_mm131_jc_mm131_fit_tmpy] = deal(zeros(interval,interval,length(x_doses)));
[result.ic90_dose_mm131_jc_fit_tmpy_10_min_val_3] = deal(zeros(interval,interval));

[p_Akt_auc_metmab_jc_mm131_fit_tmpy, ys_ic_curve_metmab_jc_mm131_fit_tmpy, p_Akt_10_min_val_metmab_jc_mm131_fit_tmpy,p_Akt_30_min_val_metmab_jc_mm131_fit_tmpy,p_Akt_60_min_val_metmab_jc_mm131_fit_tmpy] = deal(zeros(interval,interval,length(x_doses)));
[result.ic90_dose_metmab_jc_fit_tmpy_10_min_val_3] = deal(zeros(interval,interval));

[p_Akt_auc_bispec_ours_jc_mm131_fit_tmpy, ys_ic_curve_bispec_ours_jc_mm131_fit_tmpy, p_Akt_10_min_val_bispec_ours_jc_mm131_fit_tmpy,p_Akt_30_min_val_bispec_ours_jc_mm131_fit_tmpy,p_Akt_60_min_val_bispec_ours_jc_mm131_fit_tmpy] = deal(zeros(interval,interval,length(x_doses)));
[result.ic90_dose_bispec_ours_jc_fit_tmpy_10_min_val_3] = deal(zeros(interval,interval));

[p_Akt_auc_bispec_lilly_jc_mm131_fit_tmpy, ys_ic_curve_bispec_lilly_jc_mm131_fit_tmpy, p_Akt_10_min_val_bispec_lilly_jc_mm131_fit_tmpy,p_Akt_30_min_val_bispec_lilly_jc_mm131_fit_tmpy,p_Akt_60_min_val_bispec_lilly_jc_mm131_fit_tmpy] = deal(zeros(interval,interval,length(x_doses)));
[result.ic90_dose_bispec_lilly_jc_fit_tmpy_10_min_val_3] = deal(zeros(interval,interval));

[p_Akt_auc_mm151131_jc_mm131_fit_tmpy, ys_ic_curve_mm151131_jc_mm131_fit_tmpy, p_Akt_10_min_val_mm151131_jc_mm131_fit_tmpy,p_Akt_30_min_val_mm151131_jc_mm131_fit_tmpy,p_Akt_60_min_val_mm151131_jc_mm131_fit_tmpy] = deal(zeros(interval,interval,length(x_doses)));
[result.ic90_dose_mm151131_jc_fit_tmpy_10_min_val_3] = deal(zeros(interval,interval));

q2 = zeros(3,4);

tic

% Fix EGFR to 2.000.000

arSetPars('init_egfr',log10(2000000));
ar.model.qPlotYs(:) = 0;

%Range of met and epcam 1.000 -> 500.000 // -5 -> 0.4

for i = 1:interval
    arSetPars('init_met',log10(min_val_met) + (log10(max_val_met)-log10(min_val_met))/interval * i );
    for j = 1:interval
        if(mod(i,3)==0 && mod(j,3)==0)
            save tmp_inter.mat
        end
        ic90_dose(i,j) = 0;
        ic50_dose(i,j) = 0;       
        arSetPars('init_epcam',log10(min_val_ep) + (log10(max_val_ep)-log10(min_val_ep))/interval * j);
        die_ps_jc_fit_tmpy(i,j,:) = ar.p;
        
        if(single)
        %read-out for mm131  
        ar.model.qPlotXs(:) = 0;
        ar.model.qPlotXs(4) = 1;
        
        arPlot
        h = gcf;
        axesObjs = get(h,'Children');
        dataObjs = get(axesObjs,'Children');
        
        for k = 1:length(x_doses)
            p_Akt_auc_mm131_jc_mm131_fit_tmpy(i,j,k) = 0;
            ys_ic_curve_mm131_jc_mm131_fit_tmpy(i,j,k) = 0;
            pAkt_x = dataObjs{pos_akt}(k).XData;
            pAkt_y = dataObjs{pos_akt}(k).YData;
            
            if(length(pAkt_x) == 830)
                p_Akt_10_min_val_mm131_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(409);
                p_Akt_30_min_val_mm131_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(667);
                p_Akt_60_min_val_mm131_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(803);
                
                for p = 1:81
                    p_Akt_auc_mm131_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_mm131_jc_mm131_fit_tmpy(i,j,k) + 0.5*(pAkt_y((p-1)*10+1) + pAkt_y(p*10)) * (pAkt_x(p*10) - pAkt_x((p-1)*10+1));
                end
            elseif(length(pAkt_x) == 590)
                p_Akt_10_min_val_mm131_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(192);
                p_Akt_30_min_val_mm131_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(427);
                p_Akt_60_min_val_mm131_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(563);
                
                for p = 1:57
                    p_Akt_auc_mm131_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_mm131_jc_mm131_fit_tmpy(i,j,k) + 0.5*(pAkt_y((p-1)*10+1) + pAkt_y(p*10)) * (pAkt_x(p*10) - pAkt_x((p-1)*10+1));
                    % p_Akt_auc_mm131_norm_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_mm131_jc_mm131_fit_tmpy(i,j,k)/57;
                end
            else
                sprintf('New Simu_length( %i ), abording!', length(pAkt_x))
                return;
            end
            y_tmp(k) = p_Akt_10_min_val_mm131_jc_mm131_fit_tmpy(i,j,k);
        end
        
        
%         y_tmp(1) = p_Akt_10_min_val_mm131_jc_mm131_fit_tmpy(i,j,1);
%         y_tmp(2) = p_Akt_10_min_val_mm131_jc_mm131_fit_tmpy(i,j,2);
%         y_tmp(3) = p_Akt_10_min_val_mm131_jc_mm131_fit_tmpy(i,j,3);
%         y_tmp(4) = p_Akt_10_min_val_mm131_jc_mm131_fit_tmpy(i,j,4);
%         y_tmp(5) = p_Akt_10_min_val_mm131_jc_mm131_fit_tmpy(i,j,5);
%         y_tmp(6) = p_Akt_10_min_val_mm131_jc_mm131_fit_tmpy(i,j,6);
%         y_tmp(7) = p_Akt_10_min_val_mm131_jc_mm131_fit_tmpy(i,j,7);
%         y_tmp(8) = p_Akt_10_min_val_mm131_jc_mm131_fit_tmpy(i,j,8);
        
        y_tmp = sort(y_tmp);

        if(ic90_req_10_min_val < y_tmp(1))
           sprintf('ID i j %d %d (high)', i, j)
           result.ic90_dose_mm131_jc_fit_tmpy_10_min_val_3(i,j) = 4;
        elseif(ic90_req_10_min_val > y_tmp(length(x_doses)))
           sprintf('ID i j %d %d (low)', i, j)
           result.ic90_dose_mm131_jc_fit_tmpy_10_min_val_3(i,j) = -2;
        else
            for k = 1:3
                try
                    q2(k,:) = nlinfit(y_tmp,x_doses,sig_func_two,[k*0.33 k*0.33 k*0.33 k*0.33]);
                    if(isreal(q2(k,:)))
                        sprintf('ID 1 k i j %d %d %d', k, i, j)
                        result.ic90_dose_mm131_jc_fit_tmpy_10_min_val_3(i,j) = sig_func_two(q2(k,:),ic90_req_10_min_val);
                    end
                catch

                end
            end
        end

        close all
        end
        
        
        if(single)
        %read_out for metmab
        ar.model.qPlotXs(:) = 0;
        ar.model.qPlotXs(6) = 1;
        
        arPlot
        h = gcf;
        axesObjs = get(h,'Children');
        dataObjs = get(axesObjs,'Children');
        
        for k = 1:length(x_doses)
            p_Akt_auc_metmab_jc_mm131_fit_tmpy(i,j,k) = 0;
            ys_ic_curve_metmab_jc_mm131_fit_tmpy(i,j,k) = 0;
            pAkt_x = dataObjs{pos_akt}(k).XData;
            pAkt_y = dataObjs{pos_akt}(k).YData;
            if(length(pAkt_x) == 830)
                p_Akt_10_min_val_metmab_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(409);
                p_Akt_30_min_val_metmab_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(667);
                p_Akt_60_min_val_metmab_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(803);

                for p = 1:81
                    p_Akt_auc_metmab_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_metmab_jc_mm131_fit_tmpy(i,j,k) + 0.5*(pAkt_y((p-1)*10+1) + pAkt_y(p*10)) * (pAkt_x(p*10) - pAkt_x((p-1)*10+1));
                end
            elseif(length(pAkt_x) == 590)
                p_Akt_10_min_val_metmab_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(192);
                p_Akt_30_min_val_metmab_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(427);
                p_Akt_60_min_val_metmab_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(563);

                for p = 1:57
                    p_Akt_auc_metmab_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_metmab_jc_mm131_fit_tmpy(i,j,k) + 0.5*(pAkt_y((p-1)*10+1) + pAkt_y(p*10)) * (pAkt_x(p*10) - pAkt_x((p-1)*10+1));
                end
            else
                sprintf('New Simu_length %i , abording!', length(pAkt_x))
                return;
            end
            y_tmp(k) = p_Akt_10_min_val_metmab_jc_mm131_fit_tmpy(i,j,k);
        end

%         y_tmp(1) = p_Akt_10_min_val_metmab_jc_mm131_fit_tmpy(i,j,1);
%         y_tmp(2) = p_Akt_10_min_val_metmab_jc_mm131_fit_tmpy(i,j,2);
%         y_tmp(3) = p_Akt_10_min_val_metmab_jc_mm131_fit_tmpy(i,j,3);
%         y_tmp(4) = p_Akt_10_min_val_metmab_jc_mm131_fit_tmpy(i,j,4);
%         y_tmp(5) = p_Akt_10_min_val_metmab_jc_mm131_fit_tmpy(i,j,5);
%         y_tmp(6) = p_Akt_10_min_val_metmab_jc_mm131_fit_tmpy(i,j,6);
%         y_tmp(7) = p_Akt_10_min_val_metmab_jc_mm131_fit_tmpy(i,j,7);
%         y_tmp(8) = p_Akt_10_min_val_metmab_jc_mm131_fit_tmpy(i,j,8);

        y_tmp = sort(y_tmp);
        
        if(ic90_req_10_min_val < y_tmp(1))
           sprintf('ID i j %d %d (high)', i, j)
           result.ic90_dose_metmab_jc_fit_tmpy_10_min_val_3(i,j) = 4;
        elseif(ic90_req_10_min_val > y_tmp(length(x_doses)))
           sprintf('ID i j %d %d (low)', i, j)
           result.ic90_dose_metmab_jc_fit_tmpy_10_min_val_3(i,j) = -2;
        else
            for k = 1:3
                try
                    q2(k,:) = nlinfit(y_tmp,x_doses,sig_func_two,[k*0.33 k*0.33 k*0.33 k*0.33]);
                    
                    if(isreal(q2(k,:)))
                        sprintf('ID 2 k i j %d %d %d', k, i, j)
%                         if(sig_func_two(q2(k,:),ic90_req_10_min_val)<0)
%                             if(~isnumeric(ic90_dose_metmab_jc_fit_tmpy_10_min_val_3(i,j)))
%                                 ic90_dose_metmab_jc_fit_tmpy_10_min_val_3(i,j) = nan;
%                             end
%                         else
                            result.ic90_dose_metmab_jc_fit_tmpy_10_min_val_3(i,j) = sig_func_two(q2(k,:),ic90_req_10_min_val);
%                         end
                    end
                    %sprintf('%s',msgid)
                catch
                    %result.ic90_dose_mm131_jc_fit_tmpy_10_min_val(i,j) = nan;
                    %result.ic90_dose_mm131_jc_fit_tmpy_10_min_val_2(i,j) = nan;
                end
            end
        end
% 
%         %ic90
%         %ic90_req = max_val(max_val(p_Akt_auc(:,:,5)))*0.1;
%         %ic90_dose_metmab_on(i,j) = polyval(q,ic90_req);
%         ic90_dose_metmab_jc_long(i,j) = sig_func(q,ic90_req);
%        
% 
%         %ic50
%         %ic50_req = max_val(max_val(p_Akt_auc(:,:,5)))*0.5;
%         ic50_dose_metmab_kd_fit_jc_long(i,j) = polyval(q,ic50_req);

        close all
        
        end
%         %read_out for mm151
%         ar.model.qPlotXs(:) = 0;
%         ar.model.qPlotXs(6) = 1;
%         
%         arPlot
%         h = gcf;
%         axesObjs = get(h,'Children');
%         dataObjs = get(axesObjs,'Children');
%         
%         for k = 1:5
%             p_Akt_auc_mm151_jc_mm131_fit_tmpy(i,j,k) = 0;
%             ys_ic_curve_mm151_jc_mm131_fit_tmpy(i,j,k) = 0;
%             pAkt_x = dataObjs{pos_akt}(k).XData;
%             pAkt_y = dataObjs{pos_akt}(k).YData;
%             p_Akt_10_min_val_mm151_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(409);
%             p_Akt_30_min_val_mm151_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(667);
%             p_Akt_60_min_val_mm151_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(803);
%             
%             for p = 1:81
%                 p_Akt_auc_mm151_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_mm151_jc_mm131_fit_tmpy(i,j,k) + 0.5*(pAkt_y((p-1)*10+1) + pAkt_y(p*10)) * (pAkt_x(p*10) - pAkt_x((p-1)*10+1));
%                 p_Akt_auc_mm151_norm_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_mm151_jc_mm131_fit_tmpy(i,j,k)/81;
%             end
%         end
%         
%  %       for p = 1:5
%  %           ys_ic_curve(i,j,k) = p_Akt_auc(i,j,p);
%  %       end
%         
%         %         y_tmp(1) = p_Akt_auc_mm151_jc_mm131_fit(i,j,1);
% %         y_tmp(2) = p_Akt_auc_mm151_jc_mm131_fit(i,j,2);
% %         y_tmp(3) = p_Akt_auc_mm151_jc_mm131_fit(i,j,3);
% %         y_tmp(4) = p_Akt_auc_mm151_jc_mm131_fit(i,j,4);
% %         y_tmp(5) = p_Akt_auc_mm151_jc_mm131_fit(i,j,5);
% 
%         y_tmp(1) = p_Akt_10_min_val_mm151_jc_mm131_fit_tmpy(i,j,1);
%         y_tmp(2) = p_Akt_10_min_val_mm151_jc_mm131_fit_tmpy(i,j,2);
%         y_tmp(3) = p_Akt_10_min_val_mm151_jc_mm131_fit_tmpy(i,j,3);
%         y_tmp(4) = p_Akt_10_min_val_mm151_jc_mm131_fit_tmpy(i,j,4);
%         y_tmp(5) = p_Akt_10_min_val_mm151_jc_mm131_fit_tmpy(i,j,5);
% 
%         if(ic90_req_10_min_val < y_tmp(1))
%            sprintf('ID i j %d %d (high)', i, j)
%            ic90_dose_mm151_jc_fit_tmpy_10_min_val_3(i,j) = 1;
%         elseif(ic90_req_10_min_val > y_tmp(5))
%            sprintf('ID i j %d %d (low)', i, j)
%            ic90_dose_mm151_jc_fit_tmpy_10_min_val_3(i,j) = -2;
%         else
%             for k = 1:3
%                 try
%                     %q(k,:) = nlinfit(y_tmp,x_doses,sig_func,[-1+k*0.1 -1+k*0.1 -1+k*0.1]);
%                     %q2(k,:) = nlinfit(y_tmp,x_doses,sig_func_two,[-1+k*0.1 -1+k*0.1 -1+k*0.1 -1+k*0.1]);
%                     q2(k,:) = nlinfit(y_tmp,x_doses,sig_func_two,[k*0.33 k*0.33 k*0.33 k*0.33]);
%                     %[msgstr, msgid] = lastwarn;
%     %                 if(isreal(q(k,:)) && q(k,:) ~= 0)
%     %                     sprintf('%d Here %d', q(k),k)      
%     %                     result.ic90_dose_mm131_jc_fit_tmpy_10_min_val(i,j) = sig_func(q(k,:),ic90_req_10_min_val);
%     %                 end
%                     if(isreal(q2(k,:)))
%                         sprintf('ID k i j %d %d %d', k, i, j)
%                         if(sig_func_two(q2(k,:),ic90_req_10_min_val)<0)
%                             if(~isnumeric(ic90_dose_mm151_jc_fit_tmpy_10_min_val_3(i,j)))
%                                 ic90_dose_mm151_jc_fit_tmpy_10_min_val_3(i,j) = nan;
%                             end
%                         else
%                             ic90_dose_mm151_jc_fit_tmpy_10_min_val_3(i,j) = sig_func_two(q2(k,:),ic90_req_10_min_val);
%                         end
%                     end
%                     %sprintf('%s',msgid)
%                 catch
%                     %result.ic90_dose_mm131_jc_fit_tmpy_10_min_val(i,j) = nan;
%                     %result.ic90_dose_mm131_jc_fit_tmpy_10_min_val_2(i,j) = nan;
%                 end
%             end
%         end
%         
%         close all
        
        if(~single)

        %read_out for mm151&mm131
        ar.model.qPlotXs(:) = 0;
        ar.model.qPlotXs(8) = 1;
        
         arPlot
        h = gcf;
        axesObjs = get(h,'Children');
        dataObjs = get(axesObjs,'Children');
        
        for k = 1:length(x_doses)
            p_Akt_auc_mm151131_jc_mm131_fit_tmpy(i,j,k) = 0;
            ys_ic_curve_mm151131_jc_mm131_fit_tmpy(i,j,k) = 0;
            pAkt_x = dataObjs{pos_akt}(k).XData;
            pAkt_y = dataObjs{pos_akt}(k).YData;
            if(length(pAkt_x) == 830)
                p_Akt_10_min_val_mm151131_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(409);
                p_Akt_30_min_val_mm151131_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(667);
                p_Akt_60_min_val_mm151131_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(803);

                for p = 1:81
                    p_Akt_auc_mm151131_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_mm151131_jc_mm131_fit_tmpy(i,j,k) + 0.5*(pAkt_y((p-1)*10+1) + pAkt_y(p*10)) * (pAkt_x(p*10) - pAkt_x((p-1)*10+1));
                    %p_Akt_auc_mm151131_norm_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_mm151131_jc_mm131_fit_tmpy(i,j,k)/81;
                end
            elseif(length(pAkt_x) == 590)
                p_Akt_10_min_val_mm151131_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(192);
                p_Akt_30_min_val_mm151131_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(427);
                p_Akt_60_min_val_mm151131_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(563);

                for p = 1:57
                    p_Akt_auc_mm151131_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_mm151131_jc_mm131_fit_tmpy(i,j,k) + 0.5*(pAkt_y((p-1)*10+1) + pAkt_y(p*10)) * (pAkt_x(p*10) - pAkt_x((p-1)*10+1));
                    %p_Akt_auc_mm151131_norm_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_mm151131_jc_mm131_fit_tmpy(i,j,k)/81;
                end
            else
                sprintf('New Simu_length( %i ), abording!', length(pAkt_x))
                return;
            end
            y_tmp(k) = p_Akt_10_min_val_mm151131_jc_mm131_fit_tmpy(i,j,k);
        end
        
 %       for p = 1:5
 %           ys_ic_curve(i,j,k) = p_Akt_auc(i,j,p);
 %       end
        
        %         y_tmp(1) = p_Akt_auc_mm151131_jc_mm131_fit(i,j,1);
%         y_tmp(2) = p_Akt_auc_mm151131_jc_mm131_fit(i,j,2);
%         y_tmp(3) = p_Akt_auc_mm151131_jc_mm131_fit(i,j,3);
%         y_tmp(4) = p_Akt_auc_mm151131_jc_mm131_fit(i,j,4);
%         y_tmp(5) = p_Akt_auc_mm151131_jc_mm131_fit(i,j,5);

%         y_tmp(1) = p_Akt_10_min_val_mm151131_jc_mm131_fit_tmpy(i,j,1);
%         y_tmp(2) = p_Akt_10_min_val_mm151131_jc_mm131_fit_tmpy(i,j,2);
%         y_tmp(3) = p_Akt_10_min_val_mm151131_jc_mm131_fit_tmpy(i,j,3);
%         y_tmp(4) = p_Akt_10_min_val_mm151131_jc_mm131_fit_tmpy(i,j,4);
%         y_tmp(5) = p_Akt_10_min_val_mm151131_jc_mm131_fit_tmpy(i,j,5);
%         y_tmp(6) = p_Akt_10_min_val_mm151131_jc_mm131_fit_tmpy(i,j,6);
%         y_tmp(7) = p_Akt_10_min_val_mm151131_jc_mm131_fit_tmpy(i,j,7);
%         y_tmp(8) = p_Akt_10_min_val_mm151131_jc_mm131_fit_tmpy(i,j,8);

        y_tmp = sort(y_tmp);
        
        if(ic90_req_10_min_val < y_tmp(1))
           sprintf('ID i j %d %d (high)', i, j)
           result.ic90_dose_mm151131_jc_fit_tmpy_10_min_val_3(i,j) = 4;
        elseif(ic90_req_10_min_val > y_tmp(length(x_doses)))
           sprintf('ID i j %d %d (low)', i, j)
           result.ic90_dose_mm151131_jc_fit_tmpy_10_min_val_3(i,j) = -2;
        else
            for k = 1:3
                try
                    %q(k,:) = nlinfit(y_tmp,x_doses,sig_func,[-1+k*0.1 -1+k*0.1 -1+k*0.1]);
                    %q2(k,:) = nlinfit(y_tmp,x_doses,sig_func_two,[-1+k*0.1 -1+k*0.1 -1+k*0.1 -1+k*0.1]);
                    q2(k,:) = nlinfit(y_tmp,x_doses,sig_func_two,[k*0.33 k*0.33 k*0.33 k*0.33]);
                    %[msgstr, msgid] = lastwarn;
    %                 if(isreal(q(k,:)) && q(k,:) ~= 0)
    %                     sprintf('%d Here %d', q(k),k)      
    %                     result.ic90_dose_mm131_jc_fit_tmpy_10_min_val(i,j) = sig_func(q(k,:),ic90_req_10_min_val);
    %                 end
                    if(isreal(q2(k,:)))
                        sprintf('ID add k i j %d %d %d', k, i, j)
%                         if(sig_func_two(q2(k,:),ic90_req_10_min_val)<0)
%                             if(~isnumeric(ic90_dose_mm151131_jc_fit_tmpy_10_min_val_3(i,j)))
%                                 ic90_dose_mm151131_jc_fit_tmpy_10_min_val_3(i,j) = nan;
%                             end
%                         else
                            result.ic90_dose_mm151131_jc_fit_tmpy_10_min_val_3(i,j) = sig_func_two(q2(k,:),ic90_req_10_min_val);
%                         end
                    end
                    %sprintf('%s',msgid)
                catch
                    %result.ic90_dose_mm131_jc_fit_tmpy_10_min_val(i,j) = nan;
                    %result.ic90_dose_mm131_jc_fit_tmpy_10_min_val_2(i,j) = nan;
                end
            end
        end
        
        close all
        
        end
        
        if(~single)
        
        %read_out for bispec_ours
        ar.model.qPlotXs(:) = 0;
        ar.model.qPlotXs(9) = 1;
        
        arPlot
        h = gcf;
        axesObjs = get(h,'Children');
        dataObjs = get(axesObjs,'Children');
        
        for k = 1:length(x_doses)
            p_Akt_auc_bispec_ours_jc_mm131_fit_tmpy(i,j,k) = 0;
            ys_ic_curve_bispec_ours_jc_mm131_fit_tmpy(i,j,k) = 0;
            pAkt_x = dataObjs{pos_akt}(k).XData;
            pAkt_y = dataObjs{pos_akt}(k).YData;
            if(length(pAkt_x) == 830)
                p_Akt_10_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(409);
                p_Akt_30_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(667);
                p_Akt_60_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(803);

                for p = 1:81
                    p_Akt_auc_bispec_ours_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_bispec_ours_jc_mm131_fit_tmpy(i,j,k) + 0.5*(pAkt_y((p-1)*10+1) + pAkt_y(p*10)) * (pAkt_x(p*10) - pAkt_x((p-1)*10+1));
                    %p_Akt_auc_bispec_ours_norm_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_bispec_ours_jc_mm131_fit_tmpy(i,j,k)/81;
                end
            elseif(length(pAkt_x) == 590)
                p_Akt_10_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(192);
                p_Akt_30_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(427);
                p_Akt_60_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(563);

                for p = 1:57
                    p_Akt_auc_bispec_ours_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_bispec_ours_jc_mm131_fit_tmpy(i,j,k) + 0.5*(pAkt_y((p-1)*10+1) + pAkt_y(p*10)) * (pAkt_x(p*10) - pAkt_x((p-1)*10+1));
                    %p_Akt_auc_bispec_ours_norm_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_bispec_ours_jc_mm131_fit_tmpy(i,j,k)/81;
                end
            else
                sprintf('New Simu_length( %i ), abording!', length(pAkt_x))
                return;
            end
            y_tmp(k) = p_Akt_10_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,k);
        end
        
 %       for p = 1:5
 %           ys_ic_curve(i,j,k) = p_Akt_auc(i,j,p);
 %       end
        
        %         y_tmp(1) = p_Akt_auc_bispec_ours_jc_mm131_fit(i,j,1);
%         y_tmp(2) = p_Akt_auc_bispec_ours_jc_mm131_fit(i,j,2);
%         y_tmp(3) = p_Akt_auc_bispec_ours_jc_mm131_fit(i,j,3);
%         y_tmp(4) = p_Akt_auc_bispec_ours_jc_mm131_fit(i,j,4);
%         y_tmp(5) = p_Akt_auc_bispec_ours_jc_mm131_fit(i,j,5);

%         y_tmp(1) = p_Akt_10_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,1);
%         y_tmp(2) = p_Akt_10_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,2);
%         y_tmp(3) = p_Akt_10_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,3);
%         y_tmp(4) = p_Akt_10_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,4);
%         y_tmp(5) = p_Akt_10_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,5);
%         y_tmp(6) = p_Akt_10_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,6);
%         y_tmp(7) = p_Akt_10_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,7);
%         y_tmp(8) = p_Akt_10_min_val_bispec_ours_jc_mm131_fit_tmpy(i,j,8);

        y_tmp = sort(y_tmp);
        
        if(ic90_req_10_min_val < y_tmp(1))
           sprintf('ID i j %d %d (high)', i, j)
           result.ic90_dose_bispec_ours_jc_fit_tmpy_10_min_val_3(i,j) = 4;
        elseif(ic90_req_10_min_val > y_tmp(length(x_doses)))
           sprintf('ID i j %d %d (low)', i, j)
           result.ic90_dose_bispec_ours_jc_fit_tmpy_10_min_val_3(i,j) = -2;
        else
            for k = 1:3
                try
                    %q(k,:) = nlinfit(y_tmp,x_doses,sig_func,[-1+k*0.1 -1+k*0.1 -1+k*0.1]);
                    %q2(k,:) = nlinfit(y_tmp,x_doses,sig_func_two,[-1+k*0.1 -1+k*0.1 -1+k*0.1 -1+k*0.1]);
                    q2(k,:) = nlinfit(y_tmp,x_doses,sig_func_two,[k*0.33 k*0.33 k*0.33 k*0.33]);
                    %[msgstr, msgid] = lastwarn;
    %                 if(isreal(q(k,:)) && q(k,:) ~= 0)
    %                     sprintf('%d Here %d', q(k),k)      
    %                     result.ic90_dose_mm131_jc_fit_tmpy_10_min_val(i,j) = sig_func(q(k,:),ic90_req_10_min_val);
    %                 end
                    if(isreal(q2(k,:)))
                        sprintf('ID 3 k i j %d %d %d', k, i, j)
%                         if(sig_func_two(q2(k,:),ic90_req_10_min_val)<0)
%                             if(~isnumeric(ic90_dose_bispec_ours_jc_fit_tmpy_10_min_val_3(i,j)))
%                                 ic90_dose_bispec_ours_jc_fit_tmpy_10_min_val_3(i,j) = nan;
%                             end
%                         else
                            result.ic90_dose_bispec_ours_jc_fit_tmpy_10_min_val_3(i,j) = sig_func_two(q2(k,:),ic90_req_10_min_val);
%                         end
                    end
                    %sprintf('%s',msgid)
                catch
                    %result.ic90_dose_mm131_jc_fit_tmpy_10_min_val(i,j) = nan;
                    %result.ic90_dose_mm131_jc_fit_tmpy_10_min_val_2(i,j) = nan;
                end
            end
        end
        
        close all
        
        end
        
        if(~single)
            
        %read_out for bispec_lilly
        ar.model.qPlotXs(:) = 0;
        ar.model.qPlotXs(10) = 1;
        
        arPlot
        h = gcf;
        axesObjs = get(h,'Children');
        dataObjs = get(axesObjs,'Children');
        
        for k = 1:length(x_doses)
            p_Akt_auc_bispec_lilly_jc_mm131_fit_tmpy(i,j,k) = 0;
            ys_ic_curve_bispec_lilly_jc_mm131_fit_tmpy(i,j,k) = 0;
            pAkt_x = dataObjs{pos_akt}(k).XData;
            pAkt_y = dataObjs{pos_akt}(k).YData;
            
            if(length(pAkt_x) == 830)
                p_Akt_10_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(409);
                p_Akt_30_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(667);
                p_Akt_60_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(803);

                for p = 1:81
                    p_Akt_auc_bispec_lilly_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_bispec_lilly_jc_mm131_fit_tmpy(i,j,k) + 0.5*(pAkt_y((p-1)*10+1) + pAkt_y(p*10)) * (pAkt_x(p*10) - pAkt_x((p-1)*10+1));
                    %p_Akt_auc_bispec_lilly_norm_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_bispec_lilly_jc_mm131_fit_tmpy(i,j,k)/81;
                end
            elseif(length(pAkt_x) == 590)
                p_Akt_10_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(192);
                p_Akt_30_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(427);
                p_Akt_60_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,k) = pAkt_y(563);

                for p = 1:57
                    p_Akt_auc_bispec_lilly_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_bispec_lilly_jc_mm131_fit_tmpy(i,j,k) + 0.5*(pAkt_y((p-1)*10+1) + pAkt_y(p*10)) * (pAkt_x(p*10) - pAkt_x((p-1)*10+1));
                    %p_Akt_auc_bispec_lilly_norm_jc_mm131_fit_tmpy(i,j,k) = p_Akt_auc_bispec_lilly_jc_mm131_fit_tmpy(i,j,k)/81;
                end
            else
                sprintf('New Simu_length( %i ), abording!', length(pAkt_x))
                return;
            end
            y_tmp(k) = p_Akt_10_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,k);           
        end
        
 %       for p = 1:5
 %           ys_ic_curve(i,j,k) = p_Akt_auc(i,j,p);
 %       end
        
        %         y_tmp(1) = p_Akt_auc_bispec_lilly_jc_mm131_fit(i,j,1);
%         y_tmp(2) = p_Akt_auc_bispec_lilly_jc_mm131_fit(i,j,2);
%         y_tmp(3) = p_Akt_auc_bispec_lilly_jc_mm131_fit(i,j,3);
%         y_tmp(4) = p_Akt_auc_bispec_lilly_jc_mm131_fit(i,j,4);
%         y_tmp(5) = p_Akt_auc_bispec_lilly_jc_mm131_fit(i,j,5);

%         y_tmp(1) = p_Akt_10_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,1);
%         y_tmp(2) = p_Akt_10_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,2);
%         y_tmp(3) = p_Akt_10_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,3);
%         y_tmp(4) = p_Akt_10_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,4);
%         y_tmp(5) = p_Akt_10_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,5);
%         y_tmp(6) = p_Akt_10_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,6);
%         y_tmp(7) = p_Akt_10_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,7);
%         y_tmp(8) = p_Akt_10_min_val_bispec_lilly_jc_mm131_fit_tmpy(i,j,8);
        
        y_tmp = sort(y_tmp);

        if(ic90_req_10_min_val < y_tmp(1))
           sprintf('ID i j %d %d (high)', i, j)
           result.ic90_dose_bispec_lilly_jc_fit_tmpy_10_min_val_3(i,j) = 4;
        elseif(ic90_req_10_min_val > y_tmp(length(x_doses)))
           sprintf('ID i j %d %d (low)', i, j)
           result.ic90_dose_bispec_lilly_jc_fit_tmpy_10_min_val_3(i,j) = -2;
        else
            for k = 1:3
                try
                    %q(k,:) = nlinfit(y_tmp,x_doses,sig_func,[-1+k*0.1 -1+k*0.1 -1+k*0.1]);
                    %q2(k,:) = nlinfit(y_tmp,x_doses,sig_func_two,[-1+k*0.1 -1+k*0.1 -1+k*0.1 -1+k*0.1]);
                    q2(k,:) = nlinfit(y_tmp,x_doses,sig_func_two,[k*0.33 k*0.33 k*0.33 k*0.33]);
                    %[msgstr, msgid] = lastwarn;
    %                 if(isreal(q(k,:)) && q(k,:) ~= 0)
    %                     sprintf('%d Here %d', q(k),k)      
    %                     result.ic90_dose_mm131_jc_fit_tmpy_10_min_val(i,j) = sig_func(q(k,:),ic90_req_10_min_val);
    %                 end
                    if(isreal(q2(k,:)))
                        sprintf('ID 4 k i j %d %d %d', k, i, j)
%                         if(sig_func_two(q2(k,:),ic90_req_10_min_val)<0)
%                             if(~isnumeric(ic90_dose_bispec_lilly_jc_fit_tmpy_10_min_val_3(i,j)))
%                                 ic90_dose_bispec_lilly_jc_fit_tmpy_10_min_val_3(i,j) = nan;
%                             end
%                         else
                            result.ic90_dose_bispec_lilly_jc_fit_tmpy_10_min_val_3(i,j) = sig_func_two(q2(k,:),ic90_req_10_min_val);
%                         end
                    end
                    %sprintf('%s',msgid)
                catch
                    sprintf('Help!')
                    %result.ic90_dose_mm131_jc_fit_tmpy_10_min_val(i,j) = nan;
                    %result.ic90_dose_mm131_jc_fit_tmpy_10_min_val_2(i,j) = nan;
                end
            end
        end        
        close all
        end
    end
end

toc

% ar.p(4) = met;
% ar.p(2) = egfr;

save rerun_tmp.mat

r{1} = result;
r{2} = die_ps_jc_fit_tmpy;

arSetPars('init_egf',0,0,0);

end
