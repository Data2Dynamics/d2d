function f = correlationmatrixplot(A,names,p_val)
%CORRELATIONMATRIXPLOT
%
% 2015 Max Schelker

if(~exist('p_val','var'))
    p_val = 0.01;
end

scale = 10/size(A,1);
f = figure('Name','Correlation matrix');
f.Color = 'w';

colors = redwhiteblue(11);
[R,p] = corr(A');
R(isnan(R)) = 0;
p(isnan(p)) = 1;

l = 0;
for i=1:size(A,1)
    for j=1:size(A,1)
        l = l+1;
        subaxis(size(A,1),size(A,1),l, 'Spacing', 0.005, 'Padding', 0, 'Margin', 0.01);
        if i>j
            cindex = round((1+R(i,j))*5)+1;
            xl = [min(A(i,:)) max(A(i,:))];
            if(diff(xl)==0)
                xl = [0 1];
            end
            yl = [min(A(j,:)) max(A(j,:))];
            if(diff(yl)==0)
                yl = [0 1];
            end
            axis([xl yl]);
            if(p(i,j)<p_val)
                hp = patch([xl(1) xl(2) xl(2) xl(1)],[yl(1) yl(1) yl(2) yl(2)],colors(cindex,:));
                hp.FaceAlpha = 0.5;
                hp.EdgeColor = 'w';
                hold on;
            end
            scatter(A(i,:),A(j,:),'.k');
            hold off
        elseif i<j
            text(0.2,.5,sprintf('%1.2f',R(i,j)),'Units','normalized','FontSize',scale*11)
            if(p(i,j)<p_val)
                cindex = round((1+R(i,j))*5)+1;
                xl = [0 1];
                yl = [0 1];
                hp = patch([xl(1) xl(2) xl(2) xl(1)],[yl(1) yl(1) yl(2) yl(2)],colors(cindex,:));
                hp.FaceAlpha = 0.5;
                hp.EdgeColor = 'w';
            end
        else
            text(0.05,.5,linewrap(names{i},10),'Units','normalized','FontSize',scale*9,'FontWeight','bold')
        end
        axis off
    end
end
