% plot merged

function arPlotMerged(condition_link, reference_condition, saveToFile)

global ar

if(isempty(ar))
    error('please initialize by arInit')
end

if(~exist('condition_link','var') || isempty(condition_link))
    condition_link = {};
    for jm=1:length(ar.model);
        condition_link_tmp = {};
        for jc = 1:length(ar.model(jm).condition)
            condition_link_tmp{jc} = jc; %#ok<AGROW>
        end
        condition_link{jm} = condition_link_tmp; %#ok<AGROW>
    end
end
if(~exist('reference_condition','var'))
    for jm=1:length(ar.model)
        reference_condition{jm} = 1;
    end
end
if(~exist('saveToFile','var'))
    saveToFile = false;
end

standcond = 'Reference'; % Insert standard condition

xrangefac = 0.1;

figcount = 1;
for jm=1:length(ar.model)
    clabel = arSummarizeConditions(jm, standcond, reference_condition{jm});
    
    figcount = do_TC(figcount, 'x', jm, xrangefac, saveToFile, condition_link{jm}, clabel);
    figcount = do_TC(figcount, 'z', jm, xrangefac, saveToFile, condition_link{jm}, clabel);
    figcount = do_TC(figcount, 'u', jm, xrangefac, saveToFile, condition_link{jm}, clabel);
end


function clabel = arSummarizeConditions(jm, standcond, reference_condition)

global ar

% collect condition labels (TC)
clabel = cell(1,length(ar.model(jm).condition));
fprintf('m%i: %s TC\n', jm, ar.model(jm).name);
for jc = 1:length(ar.model(jm).condition)
    tmpstr = sprintf('(c%i) ',jc);
    if(jc~=reference_condition)
        jd = ar.model(jm).condition(jc).dLink(1);
        jdref = ar.model(jm).condition(reference_condition).dLink(1);
        for jdc=1:length(ar.model(jm).data(jd).condition)
            jdcref = -1;
            for jdc2=1:length(ar.model(jm).data(jdref).condition)
                if(strcmp(ar.model(jm).data(jd).condition(jdc).parameter, ...
                        ar.model(jm).data(jdref).condition(jdc2).parameter))
                    jdcref = jdc2;
                end
            end
            if(jdcref == -1 || str2double(ar.model(jm).data(jd).condition(jdc).value) ~= ...
                    str2double(ar.model(jm).data(jdref).condition(jdcref).value))
                tmpstr = [tmpstr sprintf('%s=%s ',ar.model(jm).data(jd).condition(jdc).parameter, ...
                    ar.model(jm).data(jd).condition(jdc).value)]; %#ok<AGROW>
            end
        end
    else
        tmpstr = [tmpstr standcond]; %#ok<AGROW>
    end
    clabel{jc} = strtrim(tmpstr);
    
    hasTC = false;
    for jd=ar.model(jm).condition(jc).dLink
        if(ar.model(jm).data(jd).doseresponse==0)
            hasTC = true;
        end
    end
    if(hasTC)
        fprintf('m%i: %s\n', jm, clabel{jc});
    else
        fprintf('[m%i: no TC %s]\n', jm, clabel{jc});
    end
end


% % collect condition labels (DR)
% clabelDR = cell(0,4);
% fprintf('m%i: %s DR\n', jm, ar.model(jm).name);
% for jp=1:length(ar.model(jm).plot)
%     if(ar.model(jm).plot(jp).doseresponse==1)
%         responsepar = ar.model(jm).data(ar.model(jm).plot(jp).dLink(1)).response_parameter;
%         times = [];
%         for jd=ar.model(jm).plot(jp).dLink
%             times = union(times, ar.model(jm).data(jd).tExp);
%         end
%         for jt=1:length(times)
%             if(size(clabelDR,1)>0)
%                 q = ismember(clabelDR(:,1:2), {responsepar, num2str(times(jt))});
%                 q = sum(q,2)==2;
%                 if(sum(q)==0)
%                     clabelDR{end+1,1} = responsepar; %#ok<AGROW>
%                     clabelDR{end,2} = num2str(times(jt));
%                     clabelDR{end,3} = [responsepar ' t=' num2str(times(jt))];
%                     clabelDR{end,4} = [];
%                     q = ismember(clabelDR(:,1:2), {responsepar, num2str(times(jt))});
%                     q = sum(q,2)==2;
%                 end
%             else
%                 clabelDR{end+1,1} = responsepar; %#ok<AGROW>
%                 clabelDR{end,2} = num2str(times(jt));
%                 clabelDR{end,3} = [responsepar ' t=' num2str(times(jt))];
%                 clabelDR{end,4} = [];
%                 q = 1;
%             end
%             
%             for jd=ar.model(jm).plot(jp).dLink
%                 for jc = 1:length(ar.model(jm).data(jd).condition)
%                     if(strcmp(ar.model(jm).data(jd).condition(jc).parameter, responsepar))
%                         jcondi = jc;
%                     end
%                 end
%                 dose = str2double(ar.model(jm).data(jd).condition(jcondi).value);
%                 
%                 clabelDR{q,4} = union(clabelDR{q,4}, dose);
%             end
%         end
%     end
% end
% for j=1:size(clabelDR,1)
%     fprintf('m%i: (cdr%i) %s\n', jm, j, clabelDR{j,3});
% end
% fprintf('\n');



function figcount = do_TC(figcount, fname, jm, xrangefac, saveToFile, condition_link, clabel)

global ar

all_cs = [];
for jcs = 1:length(condition_link)
    all_cs = union(all_cs, condition_link{jcs});
end

colors = lines(length(ar.model(jm).condition));

% collect data
collect_cs = cell(1,length(ar.model(jm).(fname)));
collect_cs_dr = cell(1,length(ar.model(jm).(fname)));

ctime = cell(length(ar.model(jm).condition),length(ar.model(jm).(fname)));
cdata = cell(length(ar.model(jm).condition),length(ar.model(jm).(fname)));
xdata = cell(length(ar.model(jm).condition),length(ar.model(jm).(fname)));
xrange = zeros(1,length(ar.model(jm).(fname)));

ctime_dr = cell(1,length(ar.model(jm).(fname)));
cdata_dr = cell(1,length(ar.model(jm).(fname)));

for jp=find(ar.model(jm).qPlotYs)
    fprintf('merging %s ...\n', ar.model(jm).plot(jp).name);
    for jd = ar.model(jm).plot(jp).dLink
        if(strcmp(fname,'x'))
            qPlotX = ar.model(jm).qPlotX;
        elseif(strcmp(fname,'z'))
            qPlotX = ar.model(jm).qPlotZ;
        elseif(strcmp(fname,'u'))
            qPlotX = ar.model(jm).qPlotU;
        end
        for jx=find(qPlotX)
            linkname = ['merge_linker_' fname];
            if(~isfield(ar.model(jm).data(jd), linkname) || ...
                    length(ar.model(jm).data(jd).(linkname))<jx)
                ar.model(jm).data(jd).(linkname){jx} = findinobs(ar.model(jm).data(jd).fy, ...
                    ar.model(jm).(fname){jx}, ar.model(jm).(fname));
            end
            qy = ar.model(jm).data(jd).(linkname){jx};
            
            if(sum(qy)>0)
                jc = ar.model(jm).data(jd).cLink;
                if(sum(all_cs==jc)~=0)
                    if(ar.model(jm).plot(jp).doseresponse==0)
                        collect_cs{jx} = union(collect_cs{jx}, jc);
                        if(strcmp(fname,'x'))
                            xFineSimu = ar.model(jm).condition(jc).xFineSimu(:,jx);
                        elseif(strcmp(fname,'z'))
                            xFineSimu = ar.model(jm).condition(jc).zFineSimu(:,jx);
                        elseif(strcmp(fname,'u'))
                            xFineSimu = ar.model(jm).condition(jc).uFineSimu(:,jx);
                        end
                        xrange(jx) = max([xrange(jx) range(xFineSimu)]);
                    else
                        collect_cs_dr{jx} = union(collect_cs_dr{jx}, ar.model(jm).data(jd).cLink);
                    end
                    for jy=find(qy')
                        tExp = ar.model(jm).data(jd).tExp;
                        tExpLink = ar.model(jm).data(jd).tLinkExp;
                        if(strcmp(fname,'x'))
                            xExp = ar.model(jm).condition(jc).xExpSimu(tExpLink,jx);
                        elseif(strcmp(fname,'z'))
                            xExp = ar.model(jm).condition(jc).zExpSimu(tExpLink,jx);
                        elseif(strcmp(fname,'u'))
                            xExp = ar.model(jm).condition(jc).uExpSimu(tExpLink,jx);
                        end
                        res = ar.model(jm).data(jd).res(:,jy);
                        if(ar.model(jm).plot(jp).doseresponse==0)
                            ctime{jc,jx} = [ctime{jc,jx}; tExp(:)];
                            cdata{jc,jx} = [cdata{jc,jx}; res(:)];
                            xdata{jc,jx} = [xdata{jc,jx}; xExp(:)];
                        else
                            
                        end
                    end
                end
            end
        end
    end
end




% recalc link
condition_link_new = {};
for jx = find(~cellfun(@isempty, collect_cs))
    condition_link_tmp = {};
    for j = 1:length(condition_link)
        condition_link_tmp{j} = intersect(condition_link{j}, collect_cs{jx}); %#ok<AGROW>
    end
    condition_link_tmp = condition_link_tmp(~cellfun(@isempty, condition_link_tmp));
    condition_link_new{jx} = condition_link_tmp; %#ok<AGROW>
end

% plot time course
for jx = find(~cellfun(@isempty, collect_cs))
    % time course plots
    ar.model(jm).plot_merged(jx).name = sprintf('%s_%s_TC', ar.model(jm).name, ar.model(jm).(fname){jx});
    h = myRaiseFigure(jm, jx, ar.model(jm).plot_merged(jx).name, figcount);
    
    condition_link = condition_link_new{jx};
    
    [nrows, ncols] = arNtoColsAndRows(length(condition_link));
    
    gs = [];
    for jcs = 1:length(condition_link)
        colindex = 1;
        hstmp = [];
        clabelstmp = {};
        
        g = subplot(nrows, ncols, jcs);
        gs(end+1) = g; %#ok<AGROW>
        
        for jc = condition_link{jcs}
            hold(g, 'on');
            
            tFine = ar.model(jm).condition(jc).tFine;
            qt = tFine < max(ctime{jc,jx}(:))*1.1;
            % qt = true(size(tFine));
            xFine = ar.model(jm).condition(jc).xFineSimu(qt,jx);
            if(strcmp(fname,'x'))
                xFine = ar.model(jm).condition(jc).xFineSimu(qt,jx);
            elseif(strcmp(fname,'z'))
                xFine = ar.model(jm).condition(jc).zFineSimu(qt,jx);
            elseif(strcmp(fname,'u'))
                xFine = ar.model(jm).condition(jc).uFineSimu(qt,jx);
            end
            if(sum(qt)>0)
                hstmp(end+1) = plot(g, tFine(qt), xFine, 'Color', colors(colindex,:)); %#ok<AGROW>
                
                xrangefactmp = xrange(jx)*xrangefac;
                tFineP = [tFine(qt); flipud(tFine(qt))];
                xFineP = [xFine + xrangefactmp/2; flipud(xFine - xrangefactmp/2)];
                patch(tFineP, xFineP, zeros(size(xFineP))-2, ones(size(xFineP)), ...
                    'FaceColor', colors(colindex,:)*0.1+0.9, 'EdgeColor', colors(colindex,:)*0.1+0.9)
                patch(tFineP, xFineP, zeros(size(xFineP))-1, ones(size(xFineP)), 'LineStyle', '--', ...
                    'FaceColor', 'none', 'EdgeColor', colors(colindex,:)*0.3+0.7)
                plot(g, ctime{jc,jx}, xdata{jc,jx} + xrangefactmp*cdata{jc,jx}, '*', 'Color', colors(colindex,:))
                
                clabelstmp{end+1} = strrep(clabel{jc},'_','\_'); %#ok<AGROW>
            end
            colindex = colindex + 1;
        end
        
        hold(g, 'off');
        box(g, 'on');
        
        title(g, strrep(ar.model(jm).(fname){jx},'_','\_'));
        if(jc == (nrows-1)*ncols + 1)
            xlabel(g, sprintf('%s [%s]', ar.model(jm).tUnits{3}, ar.model(jm).tUnits{2}));
        end
        if(strcmp(fname,'x'))
            xUnits = ar.model(jm).xUnits;
        elseif(strcmp(fname,'z'))
            xUnits = ar.model(jm).zUnits;
        elseif(strcmp(fname,'u'))
            xUnits = ar.model(jm).uUnits;
        end
        ylabel(g, sprintf('%s [%s]', xUnits{jx,3}, xUnits{jx,2}));
        legend(hstmp, clabelstmp);
        % legend(hstmp, clabelstmp, 'Location','Best');
        
    end
    arSpacedAxisLimits(gs);
    
    % save
    if(saveToFile)
        ar.model(jm).plot_merged(jx).savePath_FigY = mySaveFigure(h, ar.model(jm).plot_merged(jx).name);
    end
    
    figcount = figcount + 1;
end


function qy = findinobs(fy, x, xs)

xother = setdiff(xs, x);

qy = false(size(fy));
for j=1:length(fy)
    varlist = symvar(fy{j});
    qx = sum(ismember(varlist, x));
    qxother = sum(ismember(varlist, xother));
    
    if(qx==1 && qxother==0)
        qy(j) = true;
    end
end


function h = myRaiseFigure(m, jplot, figname, figcount)
global ar
openfigs = get(0,'Children');

figcolor = [1 1 1];
figdist = 0.02;

ar.model(m).plot_merged(jplot).time = now;

if(isfield(ar.model(m).plot_merged(jplot), 'fighandel_y') && ~isempty(ar.model(m).plot_merged(jplot).fighandel_y) && ...
        ar.model(m).plot_merged(jplot).fighandel_y ~= 0 && ...
        sum(ar.model(m).plot_merged(jplot).fighandel_y==openfigs)>0 && ...
        strcmp(get(ar.model(m).plot_merged(jplot).fighandel_y, 'Name'), figname))
    
    h = ar.model(m).plot_merged(jplot).fighandel_y;
    figure(h);
else
    h = figure('Name', figname, 'NumberTitle','off', ...
        'Units', 'normalized', 'Position', ...
        [0.1+((figcount-1)*figdist) 0.35-((figcount-1)*figdist) 0.3 0.45]);
    set(h,'Color', figcolor);
    ar.model(m).plot_merged(jplot).fighandel_y = h;
end
clf;

function savePath = mySaveFigure(h, name)

savePath = [arSave '/FiguresMerged/Y'];

if(~exist(savePath, 'dir'))
	mkdir(savePath)
end

savePath = mypath([savePath '/' name]);

saveas(h, savePath, 'fig');
print('-depsc', savePath);
% print('-dpng', savePath);
if(ispc)
    print('-dpdf', savePath);
elseif(ismac)
    system(['/usr/local/bin/ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);
else
    system(['export LD_LIBRARY_PATH=""; ps2pdf  -dEPSCrop ' savePath '.eps '  savePath '.pdf']);
end
% plot2svg([savePath '.svg'], h);

function str = mypath(str)
str = strrep(str, ' ', '\ ');
str = strrep(str, '(', '\(');
str = strrep(str, ')', '\)');

