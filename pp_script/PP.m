clear all
close all
clc

% upload di dati che ci interessano
n_rotors = 2;
n_blades = 3;
n_elements = 10;
n_panels = 2;
n_new_particles = n_blades * (n_elements + 1);
R = [0.66 0.66];
Vtip = [ 82.938, 82.938];
t_f = 0.3;
dt = 0.5e-3;


blue = [0, 0, 1];
white = [255, 255, 255]/255;
%colors_p = [linspace(red(1),white(1),(n_blades+4))', linspace(red(2),white(2),(n_blades+4))', linspace(red(3),white(3),(n_blades+4))'];

colors_p = [1  0  0;
            1  0  1;
            0  0  1;
            0  0  0];

if(~exist('./plot'))
    mkdir './plot'
end

if(~exist('./plot/frame/'))
    mkdir './plot/frame/'
end

if(~exist('./plot/video'))
    mkdir './plot/video'
end

n_timestep = t_f / dt;

%% Acquiring data
x_wake_panel_plot = [];
y_wake_panel_plot = [];
z_wake_panel_plot = [];
particles_coord = [];

for i = 2:2:n_timestep


        blade_coord = importdata(strcat('./solution/blade_coord/coord', num2str(i), '.dat'));
        [x_blade_plot, y_blade_plot, z_blade_plot] = data2plot(blade_coord(:, 1), blade_coord(:, 2), blade_coord(:, 3), 4);

        min_x = min(min(x_blade_plot));
        min_y = min(min(y_blade_plot));
        min_z = min(min(z_blade_plot));
        max_x = max(max(x_blade_plot));
        max_y = max(max(y_blade_plot));
        max_z = max(max(z_blade_plot));

    if i < n_panels + 3 && i > 1

        wake_panel_coord = importdata(strcat('./solution/wake_panel_coord/coord', num2str(i), '.dat'));
        [x_wake_panel_plot, y_wake_panel_plot, z_wake_panel_plot] = data2plot(wake_panel_coord(:, 1), wake_panel_coord(:, 2), wake_panel_coord(:, 3), 4);

    elseif i >= n_panels + 3

        wake_panel_coord = importdata(strcat('./solution/wake_panel_coord/coord', num2str(i), '.dat'));
        [x_wake_panel_plot, y_wake_panel_plot, z_wake_panel_plot] = data2plot(wake_panel_coord(:, 1), wake_panel_coord(:, 2), wake_panel_coord(:, 3), 4);
        particles_coord = importdata(strcat('./solution/particle_coord/coord', num2str(i), '.dat'));

    end

    if(~isempty(x_wake_panel_plot) && isempty(particles_coord))
            min_x = min([min(x_blade_plot), min(x_wake_panel_plot), min_x]);
            min_y = min([min(y_blade_plot), min(y_wake_panel_plot), min_y]);
            min_z = min([min(z_blade_plot), min(z_wake_panel_plot), min_z]);
            max_x = max([max(x_blade_plot), max(x_wake_panel_plot), max_x]);
            max_y = max([max(y_blade_plot), max(y_wake_panel_plot), max_y]);
            max_z = max([max(z_blade_plot), max(z_wake_panel_plot), max_z]);
        elseif(~isempty(particles_coord))
            min_x = min([min(x_blade_plot), min(x_wake_panel_plot), min(particles_coord(:, 1)), min_x]);
            min_y = min([min(y_blade_plot), min(y_wake_panel_plot), min(particles_coord(:, 2)), min_y]);
            min_z = min([min(z_blade_plot), min(z_wake_panel_plot), min(particles_coord(:, 3)), min_z]);
            max_x = max([max(x_blade_plot), max(x_wake_panel_plot), max(particles_coord(:, 1)), max_x]);
            max_y = max([max(y_blade_plot), max(y_wake_panel_plot), max(particles_coord(:, 2)), max_y]);
            max_z = max([max(z_blade_plot), max(z_wake_panel_plot), max(particles_coord(:, 3)), max_z]);
        end

end

%%

ct = importdata('solution/CT.dat');
Omega =  Vtip(1) / R(1);
time_for_a_round = 2 * pi / Omega;
iter_temp_a_round = floor(time_for_a_round / dt);
figure(1)
for i = 1:n_rotors
    ct_vec(i, :) = ct((i-1)*n_timestep + 1: i * n_timestep);
end
ct_vec(2, :) =  -ct_vec(2, :);
plot(dt:dt:t_f, ct_vec(1, :), 'b', dt:dt:t_f, ct_vec(2, :), 'r')
title('Andamento $C_T$ nel tempo della simulazione', 'Interpreter', 'latex')
legend('Rotore Superiore', 'Rotore Inferiore', 'Interpreter', 'latex')

xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('$C_T$', 'Interpreter', 'latex')

ct_round = [];
figure(2)
for j = 1:size(ct_vec, 1)
    for i = 1:floor((size(ct_vec, 2)/iter_temp_a_round))
        ct_round(j, i) = mean( (ct_vec(j, (1 + (i-1) * iter_temp_a_round) : i * iter_temp_a_round))); 
    end
    plot(1:size(ct_round, 2), ct_round(j, :))
hold on
end

title('Andamento valor medio del ct con il numero di rotazioni')
legend('Rotore Superiore', 'Rotore Inferiore', 'Interpreter', 'latex')
xlabel('Rotazione')
ylabel('ct mediato')

%%
cp = importdata('solution/CP.dat');
Omega =  Vtip(1) / R(1);
time_for_a_round = 2 * pi / Omega;
iter_temp_a_round = floor(time_for_a_round / dt);
figure(3)
for i = 1:n_rotors
    cp_vec(i, :) = cp((i-1)*n_timestep + 1: i * n_timestep);
end
plot(dt:dt:t_f, cp_vec(1, :), 'b', dt:dt:t_f, cp_vec(2, :), 'r')
title('Andamento $C_P$ nel tempo della simulazione', 'Interpreter', 'latex')
legend('Rotore Superiore', 'Rotore Inferiore', 'Interpreter', 'latex')

xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('$C_P$', 'Interpreter', 'latex')

cp_round = [];
figure(4)
for j = 1:size(cp_vec, 1)
    for i = 1:floor((size(cp_vec, 2)/iter_temp_a_round))
        cp_round(j, i) = mean( (cp_vec(j, (1 + (i-1) * iter_temp_a_round) : i * iter_temp_a_round))); 
    end
    plot(1:size(cp_round, 2), cp_round(j, :))
hold on
end

title('Andamento valor medio del cp con il numero di rotazioni')
legend('Rotore Superiore', 'Rotore Inferiore', 'Interpreter', 'latex')
xlabel('Rotazione')
ylabel('cp mediato')


%% Plotting

number = 50;

x_wake_panel_plot = [];
particles_coord = [];
t_vec = dt:dt:t_f;

for i = 2:2:n_timestep

    t_plot_1 = strcat("time = ", num2str(t_vec(i)), "s");

    blade_coord = importdata(strcat('./solution/blade_coord/coord', num2str(i), '.dat'));
    x_blade_plot = [];
    y_blade_plot = [];
    z_blade_plot = [];
    for j = 1:n_rotors
        coord = blade_coord((j-1) * n_blades * n_elements * 4 + 1: j *n_blades * n_elements * 4 , :);
        [x_blade_plot(:, j), y_blade_plot(:, j), z_blade_plot(:, j)] = data2plot(coord(:, 1), coord(:, 2), coord(:, 3), 4);
    end

    x_wake_panel_plot = [];
    y_wake_panel_plot = [];
    z_wake_panel_plot = [];

    if i <= n_panels + 1 && i > 1

        wake_panel_coord = importdata(strcat('./solution/wake_panel_coord/coord', num2str(i), '.dat'));
        for j = 1:n_rotors
            coord = wake_panel_coord((j-1) * n_blades * n_elements * 4 * (i-1) + 1: j *n_blades * 4 * n_elements* (i-1) , :);
            [x_wake_panel_plot(:, j), y_wake_panel_plot(:, j), z_wake_panel_plot(:, j)] = data2plot(coord(:, 1), coord(:, 2), coord(:, 3), 4);
        end
        
    elseif i >= n_panels + 2

        wake_panel_coord = importdata(strcat('./solution/wake_panel_coord/coord', num2str(i), '.dat'));
        for j = 1:n_rotors
            coord = wake_panel_coord((j-1) * n_blades * 4 * n_elements * (n_panels+1) + 1: j *n_blades * 4 * n_elements* (n_panels+1) , :);
            [x_wake_panel_plot(:, j), y_wake_panel_plot(:, j), z_wake_panel_plot(:, j)] = data2plot(coord(:, 1), coord(:, 2), coord(:, 3), 4);
        end
        particles_coord = importdata(strcat('./solution/particle_coord/coord', num2str(i), '.dat'));

    end

        h = figure(1);
        plot3(x_blade_plot, y_blade_plot, z_blade_plot, 'g');
        if(~isempty(x_wake_panel_plot) && isempty(particles_coord))
            hold on
            plot3(x_wake_panel_plot, y_wake_panel_plot, z_wake_panel_plot, 'k');
            hold off

        elseif(~isempty(particles_coord))
            hold on
            plot3(x_wake_panel_plot, y_wake_panel_plot, z_wake_panel_plot, 'k');
            % Tutte  le particelle di un colore
%             scatter3(particles_coord(:, 1), particles_coord(:, 2), particles_coord(:, 3), 4, 'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75]);
            % Colori diversi a seconda del rotore
                coord1 = particles_coord(1:n_blades * (n_elements+1) * (i-3), :);
                coord2 = particles_coord(n_blades * (n_elements+1) * (i-3) + 1:2*n_blades * (n_elements+1) * (i-3), :);
                scatter3(coord1(:, 1), coord1(:, 2), coord1(:, 3), 4, 'MarkerEdgeColor','b',...
           'MarkerFaceColor',[0 .75 .75]);
                scatter3(coord2(:, 1), coord2(:, 2), coord2(:, 3), 4, 'MarkerEdgeColor','r',...
           'MarkerFaceColor',[0 .75 .75]);
            hold off

        end

        axis([min_x max_x   min_y max_y   min_z max_z+0.05]);
        daspect([1 1 1])
        text(max_x-0.5, 1, max_z , t_plot_1)
        xlabel('x')
        ylabel('y')
        zlabel('z')
%          view([90 0])
%         view([90 -90])
         view([0 0])
%         view([-0.8355, -0.7, 1])
         saveas(h,sprintf('./plot/frame/frame%d.png',i));

end


video = VideoWriter('./plot/video/visual_wake.avi');
video.FrameRate=15;
open(video);
for j = 2:2 : n_timestep

    if(exist(strcat('./plot/frame/frame', num2str(j), '.png')))

      filename = strcat('./plot/frame/frame', num2str(j), '.png');
      thisimage = imread(filename);
      writeVideo(video, thisimage);
      delete(filename)
    end

end

close(video);


