function [x_plot, y_plot, z_plot] = data2plot(x, y, z, n_point)
    
    x_plot = reshape(x, [n_point, length(x) / n_point]);
    x_plot(n_point+1, :) = x_plot(1, :);
    x_plot(n_point+2,:) = nan;
    x_plot = x_plot(:);
    y_plot = reshape(y, [n_point, length(y) / n_point]);
    y_plot(n_point+1, :) = y_plot(1, :);
    y_plot(n_point+2,:) = nan;
    y_plot = y_plot(:);
    z_plot = reshape(z, [n_point, length(z) / n_point]);
    z_plot(n_point+1, :) = z_plot(1, :);
    z_plot(n_point+2,:) = nan;
    z_plot = z_plot(:);