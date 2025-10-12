function draw_arrow(ax, start, ending, color, label)
    % Draw an arrow from start to end
    quiver3(ax, start(1), start(2), start(3), ...
        ending(1)-start(1), ending(2)-start(2), ending(3)-start(3), ...
            'Color', color, 'MaxHeadSize', 0.5, 'LineWidth', 1.5);
    % Draw a text label
    text(ax, ending(1)+0.05, ending(2), ending(3), label, 'Color', color, 'FontSize', 25);
end