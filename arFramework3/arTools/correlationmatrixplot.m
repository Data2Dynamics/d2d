function f = correlationmatrixplot(A,names)
%CORRELATIONMATRIXPLOT
%
% 2015 Max Schelker

scale = 10/size(A,1);
f = figure('Name','Correlation matrix');
l = 0;
colors = redwhiteblue(11);
for i=1:size(A,1)
    for j=1:size(A,1)
        l = l+1;
        subaxis(size(A,1),size(A,1),l, 'Spacing', 0.005, 'Padding', 0, 'Margin', 0.01);
        if i>j
            R = corr(A(i,:)',A(j,:)');
            cindex = round((1+R)*5)+1;
            scatter(A(i,:),A(j,:),'k.');
            xl = xlim;
            yl = ylim;
            cla;
            p = patch([xl(1) xl(2) xl(2) xl(1)],[yl(1) yl(1) yl(2) yl(2)],colors(cindex,:));
            p.FaceAlpha = 0.5;
            p.EdgeColor = 'w';
            hold on;
            scatter(A(i,:),A(j,:),'.k');
            hold off
            axis off
        elseif i<j
            [R,p] = corr(A(i,:)',A(j,:)');
            text(0.1,.5,sprintf('%1.2f\n%1.1g',R,p),'Units','normalized','FontSize',scale*11)
            cindex = round((1+R)*5)+1;
            xl = [0 1];
            yl = [0 1];
            p = patch([xl(1) xl(2) xl(2) xl(1)],[yl(1) yl(1) yl(2) yl(2)],colors(cindex,:));
            p.FaceAlpha = 0.5;
            p.EdgeColor = 'w';
            axis off
        else
            text(0.05,.5,linewrap(names{i},10),'Units','normalized','FontSize',scale*9,'FontWeight','bold')
            axis off
        end
    end
end

