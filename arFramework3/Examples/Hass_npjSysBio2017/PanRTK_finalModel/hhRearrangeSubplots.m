%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                           %%%
%%%   Setup                                   %%%
%%%   Helge Hass, 2015                        %%%
%%%   email: helge.hass@fdm.uni-freiburg.de   %%%
%%%                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Script to re-arrange subplots, delete single plots and save result
%
%Before using Rearrange, look at ar.model(m).plot.dLink to get data links
%for the plots.
%Then, for the (one) plot of arPlotter, look at the first dLink in
%ar.model(m).data(d).yNames to get a list of names. After the first name,
%the legend is plotted!

function hhRearrangeSubplots(order,saveToFile,trunc,name_append,figh)
% Adapted from: http://www.mathworks.com/matlabcentral/newsreader/view_thread/136971

if(~exist('figh','var') | figh==0)
	figh = get(0,'Children');
end
if(~exist('saveToFile','var'))
	saveToFile = false;
end
if(~exist('trunc','var'))
	trunc=false;
end
if(~exist('name_append','var'))
	name_append=[];
end

myfigh = figh'
    h = get(myfigh, 'Children');
    save_name = myfigh.Name;
    save_name = strrep(save_name,'Y: ','');
    save_name = strrep(save_name,'X: ','');
    save_name = [save_name '_' name_append];
    h = h(setdiff(1:length(h),strcmp('legend',get(h,'Tag')))); % Do not count legends as subplots
    if(~exist('order','var'))
        order = 1:length(h);
    elseif (length(order) < length(h) & ~trunc)
        order(end+1:length(h)) = (length(order)+1):length(h);
    end
    
    fig_new = figure;
    
    xfac = .8;
    yfac = 1.;
    fontsize = 8;
    
    rowstocols = .42;
if(trunc)
        r = round(length(order)^rowstocols);
    	c = ceil(length(order) / r);
        g = h;
        
        order=length(h)-order+1;
        
        inv_order=1:length(h);
        inv_order(order)=0;
        inv_order=sort(inv_order);
else
    	r = ceil(length(h)^rowstocols);
    	c = ceil(length(h) / r);
    	g = h;
end
    for i=1:length(h)
       figure(fig_new);
        if(i>length(order) & trunc)
            
            continue
        end
       
       g(order(i)) = subplot(r, c, i);
      
    end
    
    
    for i=1:length(h)
        figure(fig_new);
        if(i>length(order) & trunc)
           delete(h(inv_order(i)));
           continue
        end
        
        oldpos=get(g(order(i)), 'pos');
        if(i>c)
            set(h(order(i)), 'pos', [oldpos(1) 0.85*oldpos(2) oldpos(3) oldpos(4)]); % Repositioning
        else
            set(h(order(i)), 'pos', [oldpos(1) oldpos(2) oldpos(3) oldpos(4)]); % Repositioning
        end
           
        
    end
    delete(fig_new);

if saveToFile
    for myfigh = figh'
        arSaveFigure(myfigh, save_name,'/Figures/Y')
    end
end
