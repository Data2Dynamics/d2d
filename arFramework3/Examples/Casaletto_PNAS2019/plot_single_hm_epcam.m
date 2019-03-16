function out = plot_single_hm_epcam(M,min_val,max_val,x_dim,y_dim,mm)
    if(~exist('x_dim','var'))
        x_dim = length(M);
    end
    if(~exist('y_dim','var'))
        y_dim = length(M);
    end
    if(~exist('min_val','var'))
        min_val = 6.3e2;
    end
    if(~exist('max_val','var'))
        max_val = 5e5;
    end
    if(~exist('mm','var'))
        %print('hello');
    end
    
    if(~exist('mm','var'))
%         M_t = transpose(M(1:x_dim,1:y_dim));
        M_t = M;
        fig_hm = figure('units','normalized','position',[.1 .1 0.6 0.8]);
        map = colormap('jet');
        map([1 length(map)],1:3) = 0.8;
        colormap(map);
        %Change out of scale values to grey
        imagesc((M_t));
        colorbar;
        ylabel('EpCam Receptor #')
        xlabel('Met Receptor #')
        %axis([6e2 5e5 6e2 5e5])
        NumTicks = 10;
        L = get(gca,'XLim');
        set(gca,'XTick',linspace(L(1),L(2),NumTicks));
        xt=get(gca,'xticklabel');
        %Dynamically adjust scale
        
        %set(gca,'xticklabel',[6.3e2 1.3e3 2.8e3 5.8e3 1.2e4 2.6e4 5.4e4 1.1e5 2.4e5 5e5]);
     
        set(gca,'xticklabel',[min_val 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*2) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*3) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*4) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*5) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*6) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*7) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*8) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*9) max_val]);
        
        %
        
        NumTicks = 10;
        L = get(gca,'YLim');
        set(gca,'YTick',linspace(L(1),L(2),NumTicks));
        yt=get(gca,'yticklabel');
        %set(gca,'yticklabel',[1.3e3 2.8e3 5.8e3 1.2e4 2.6e4 5.4e4 1.1e5 2.4e5 5e5 1e6]);
        %set(gca,'yticklabel',[6.3e2 1.3e3 2.8e3 5.8e3 1.2e4 2.6e4 5.4e4 1.1e5 2.4e5 5e5]);
        set(gca,'yticklabel',[min_val 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*2) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*3) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*4) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*5) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*6) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*7) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*8) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*9) max_val]);
        
        
        %set(gca,'yticklabel',[1.1e3 2.1e3 3.8e3 7e3 1.3e4 2.3e4 4.2e4 7.7e4 1.4e5 2.5e5]);
        caxis([-2,2])
        %caxis([0,0.12])
        set(gca,'ydir','normal')
        %caxis([min(min(M_t)),max(max(M_t))])
       % saveas(fig_hm,sprintf('heatmaps_paper2/p_Akt_ic90_%s.png',ID));
       % saveas(fig_hm,sprintf('heatmaps_paper2/p_Akt_ic90_%s.fig',ID));
    else
        sub_coords(60,60) = 0;
%         M_t = transpose(M(1:x_dim,1:y_dim));
%         M2_t = transpose(mm(1:x_dim,1:y_dim));

        M_t = M;
        M2_t = mm;

        fig_hm_diff = figure('units','normalized','position',[.1 .1 0.6 0.8]);
        map = colormap('jet');
        map([1 length(map)],1:3) = 0.8;
        colormap(map);
        for i = 1:60
            for j = 1:60
                if((M_t(i,j) == -2 && M2_t(i,j) ~= -2 ))
                    sub_coords(i,j) = -2;
                elseif (M_t(i,j) ~= 4 && M2_t(i,j) == 4)
                    sub_coords(i,j) = 4;
                else
                    sub_coords(i,j) = 0;
                end
            end
        end
        looky = (M2_t-M_t);
        for i = 1:60
            for j = 1:60
                if(sub_coords(i,j) == -2)
                    for k = 1:60
                        if(sub_coords(i,j+k) ~= -2)
                            looky(i,j) = looky(i,j+k);
                            break
                        end
                    end
                elseif(sub_coords(i,j) == 4)
                    for k = 1:60
                        if(sub_coords(i,j-k) ~= 4)
                            looky(i,j) = looky(i,j-k);
                            break
                        end
                    end
                end
            end
        end
        for i = 1:60
            for j = 1:60
                if(looky(i,j) == 0)
                    looky(i,j) = -2;
                end
            end
        end
        imagesc(looky);
        colorbar;
        ylabel('EpCam Receptor #')
        xlabel('Met Receptor #')
        %axis([6e2 5e5 6e2 5e5])
        NumTicks = 10;
        L = get(gca,'XLim');
        set(gca,'XTick',linspace(L(1),L(2),NumTicks));
        xt=get(gca,'xticklabel');
        %set(gca,'xticklabel',[1.3e3 2.8e3 5.8e3 1.2e4 2.6e4 5.4e4 1.1e5 2.4e5 5e5 1e6]);
        %set(gca,'xticklabel',[1.1e3 2.1e3 3.8e3 7e3 1.3e4 2.3e4 4.2e4 7.7e4 1.4e5 2.5e5]);
        %set(gca,'xticklabel',[6.3e2 1.3e3 2.8e3 5.8e3 1.2e4 2.6e4 5.4e4 1.1e5 2.4e5 5e5]);
        set(gca,'xticklabel',[min_val 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*2) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*3) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*4) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*5) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*6) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*7) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*8) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*9) max_val]);
        
        %set(gca,'xticklabel',[1.6e3 4.2e3 1.1e4 2.9e4 7.4e4 1.9e5]);
        NumTicks = 10;
        L = get(gca,'YLim');
        set(gca,'YTick',linspace(L(1),L(2),NumTicks));
        yt=get(gca,'yticklabel');
        %set(gca,'yticklabel',[1.3e3 2.8e3 5.8e3 1.2e4 2.6e4 5.4e4 1.1e5 2.4e5 5e5 1e6]);
        %set(gca,'yticklabel',[6.3e2 1.3e3 2.8e3 5.8e3 1.2e4 2.6e4 5.4e4 1.1e5 2.4e5 5e5]);
        
        set(gca,'yticklabel',[min_val 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*2) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*3) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*4) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*5) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*6) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*7) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*8) 10^(log10(min_val)+(log10(max_val)-log10(min_val))/10*9) max_val]);
        %set(gca,'yticklabel',[1.1e3 2.1e3 3.8e3 7e3 1.3e4 2.3e4 4.2e4 7.7e4 1.4e5 2.5e5]);
        set(gca,'ydir','normal')
        %caxis([-6,6])
        caxis([-2,2])
        out = looky;
        %saveas(fig_hm_diff,sprintf('heatmaps_paper2/p_Akt_ic90_diff_%s.png',ID));
        %saveas(fig_hm_diff,sprintf('heatmaps_paper2/p_Akt_ic90_diff_%s.fig',ID));
    end
end