function qBar( labels, y, yExp, ystdExpSimu, colors, baseValue )

    global errorBarMode;
    switch errorBarMode
        case 1
            ax = axes;
            j = bar( 1 : numel(y(:,1)), y ); hold on;
            mini = min([min(min(y)),min(min(yExp - ystdExpSimu))]);
            maxi = max([max(max(y)),max(max(yExp + ystdExpSimu))]);

            if isempty( baseValue )
                set( j(1), 'BaseValue', mini );
            else
                set( j(1), 'BaseValue', baseValue );
            end

            set(j(1), 'FaceColor', [1, 1, 1], 'EdgeColor', colors{1} );
            set(j(2), 'FaceColor', [1, 1, 1], 'EdgeColor', colors{2} );
            set(j(3), 'FaceColor', [1, 1, 1], 'EdgeColor', colors{3} );

            dx = .225;
            for b = 1 : 3
                errorbar( [ 1+(b-2)*dx; 2+(b-2)*dx; 3+(b-2)*dx ], yExp(:,b).', abs((yExp(:,b) - ystdExpSimu(:,b))-yExp(:,b)), abs((yExp(:,b) + ystdExpSimu(:,b))-yExp(:,b)), 'o', 'Color', colors{b}, 'LineWidth', .5, 'MarkerSize', 3 );
            end

            set(ax, 'xTickLabel', labels);
            set(gca,'TickDir','out');
            xlim([.5, 3.5]);
            ylim([mini, maxi+0.1*(maxi-mini)]);

            set(ax,'YGrid','on');
            box off;
        case 2
            ax = axes;
            
            % Plot the bars anyway, since they set up the axes the way we
            % want (bit of a hack?)
            j = bar( 1 : numel(y(:,1)), y ); hold on;
            mini = min([min(min(y)),min(min(yExp - ystdExpSimu))]);
            maxi = max([max(max(y)),max(max(yExp + ystdExpSimu))]);

            if isempty( baseValue )
                set( j(1), 'BaseValue', mini );
            else
                set( j(1), 'BaseValue', baseValue );
            end

            set(j, 'Visible', 'off');
            
            dx = .225;
            dx2 = .05;
            for b = 1 : 3
                errorbar( [ 1+(b-2)*dx; 2+(b-2)*dx; 3+(b-2)*dx ]-dx2, y(:,b).', abs((y(:,b) - ystdExpSimu(:,b))-y(:,b)), abs((y(:,b) + ystdExpSimu(:,b))-y(:,b)), 's', 'Color', colors{b}, 'LineWidth', .5, 'MarkerSize', 3 );
                plot( [ 1+(b-2)*dx; 2+(b-2)*dx; 3+(b-2)*dx ]+dx2, yExp(:,b).', 'o', 'Color', colors{b}, 'LineWidth', .5, 'MarkerSize', 3 );
            end

            set(ax, 'xTickLabel', labels);
            set(gca,'TickDir','out');
            xlim([.5, 3.5]);
            ylim([mini, maxi+0.1*(maxi-mini)]);

            set(ax,'YGrid','on');
            box off;

    end
end