clear all
close all
clc

fileID = fopen('vortex_particle1.in');

% upload di dati che ci interessano
n_rotors = cell2mat(textscan(fileID, '%*s %*s %d',1, 'HeaderLines', 2));
n_blades = cell2mat(textscan(fileID, '%*s %*s %d',1, 'HeaderLines', 1));
n_elements = cell2mat(textscan(fileID, '%*s %*s %d',1, 'HeaderLines', 1));
n_panels = cell2mat(textscan(fileID, '%*s %*s %d',1, 'HeaderLines', 1));
n_new_particles = n_blades * (n_elements + 1);
fileID = fopen('vortex_particle2.in');
R = cell2mat(textscan(fileID, '%*s %*s %f',1, 'HeaderLines', 3));
Vtip = cell2mat(textscan(fileID, '%*s %*s %f',1, 'HeaderLines', 4));
t_f = cell2mat(textscan(fileID, '%*s %*s %f',1, 'HeaderLines', 9));
dt = cell2mat(textscan(fileID, '%*s %*s %f',1, 'HeaderLines', 1));

% SEGUENTI 2 RIGHE DA USARE NEL CASO SI MANOMETTA IL FILE DI INPUT
% AGGIUNGENDO RIGHE A CASO
% R = 0.8;
% Vtip = 100,53;
t_f = 0.3;
dt = 0.5e-3;
% t_f = 0.285;

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
Omega = 100.53 / 0.8;
time_for_a_round = 2 * pi / Omega;
iter_temp_a_round = floor(time_for_a_round / dt);
figure(3)
for i = 1:n_rotors
    ct_vec(i, :) = ct((i-1)*n_timestep + 1: i * n_timestep);
end
% ct_vec(2, :) =  -ct_vec(2, :);
plot(dt:dt:t_f, ct_vec(1, :), 'b')%, ...
%     dt:dt:t_f, ct_vec(2, :), 'r')
title('Andamento $C_T$ nel tempo della simulazione', 'Interpreter', 'latex')
% legend('Rotore Superiore', 'Rotore Inferiore', 'Interpreter', 'latex')

xlabel('Time [s]', 'Interpreter', 'latex')
ylabel('$C_T$', 'Interpreter', 'latex')

ct_round = [];
figure(4)
for j = 1:size(ct_vec, 1)
    for i = 1:floor((size(ct_vec, 2)/iter_temp_a_round))
        ct_round(j, i) = mean( (ct_vec(j, (1 + (i-1) * iter_temp_a_round) : i * iter_temp_a_round))); 
    end
    plot(1:size(ct_round, 2), ct_round(j, :))
hold on
end

title('Andamento valor medio del ct con il numero di rotazioni')
% legend('Rotore Superiore', 'Rotore Inferiore', 'Interpreter', 'latex')
xlabel('Rotazione')
ylabel('ct mediato')

%% Altri dati

comp_time_vi = importdata('./solution/comp_time_vi.dat');
comp_time_gamma = importdata('./solution/comp_time_gamma.dat');
num_particles = importdata('./solution/num_particles.dat');

figure(5)
plot(1:length(comp_time_vi), comp_time_vi ./ max(comp_time_vi), 'b', ...
     1:length(comp_time_vi), (num_particles.*num_particles) ./ max(num_particles.*num_particles), 'r')
legend('Computational time', 'N^2')
xlabel('Iterazione temporale')
ylabel('tempo di calcolo')
title('Andamento Tempo di calcolo per la vi con il numero di particelle')

elapsed_time = sum(comp_time_vi) / 8 / 60


%% Plotting

number = 50;

% colors_particle_height = ([linspace(blue(1),white(1),number)', linspace(blue(2),white(2),number)', linspace(blue(3),white(3),number)']);
% height_discretization = linspace(min_z, max_z, number);
% toll = abs(max_z - min_z) / number;

x_wake_panel_plot = [];
particles_coord = [];
t_vec = dt:dt:t_f;

% video = VideoWriter('./plot/video/visual15.avi');
% video.FrameRate=15;
% open(video);

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
%         wake_panel_coord_old = wake_panel_coord;
%         [x_wake_panel_plot_old, y_wake_panel_plot_old, z_wake_panel_plot_old] = data2plot(wake_panel_coord_old(:, 1), wake_panel_coord_old(:, 2), wake_panel_coord_old(:, 3), 4);
%         blade_coord_old = blade_coord;
%         [x_blade_plot_old, y_blade_plot_old, z_blade_plot_old] = data2plot(blade_coord_old(:, 1), blade_coord_old(:, 2), blade_coord_old(:, 3), 4);
    elseif i >= n_panels + 2

        wake_panel_coord = importdata(strcat('./solution/wake_panel_coord/coord', num2str(i), '.dat'));
        for j = 1:n_rotors
            coord = wake_panel_coord((j-1) * n_blades * 4 * n_elements * (n_panels+1) + 1: j *n_blades * 4 * n_elements* (n_panels+1) , :);
            [x_wake_panel_plot(:, j), y_wake_panel_plot(:, j), z_wake_panel_plot(:, j)] = data2plot(coord(:, 1), coord(:, 2), coord(:, 3), 4);
        end
        particles_coord = importdata(strcat('./solution/particle_coord/coord', num2str(i), '.dat'));
%         particles_coord2 = importdata(strcat('./solution/particle_coords/coords', num2str(i), '.dat'));
%         h = figure(2);
%         plot3(x_blade_plot, y_blade_plot, z_blade_plot, 'b');
%             hold on
%             vec = [];
%             for j = 1:n_blades * n_elements
%                 ((7+ 12 * (j-1)):(12+ 12 * (j-1)))';
%                 vec = [vec; ((7+ 18 * (i-1)):(13+ 18 * (i-1)))'];
%                 vec = ((7+ 12 * (j-1)):(12+ 12 * (j-1)))';
%                 vec_old = ((1+ 12 * (j-1)):(6+ 12 * (j-1)))';
%                 plot3(x_wake_panel_plot(vec, :), y_wake_panel_plot(vec, :), z_wake_panel_plot(vec, :), 'r');
%                 plot3(x_wake_panel_plot_old(vec_old, :), y_wake_panel_plot(vec_old, :), z_wake_panel_plot(vec_old, :), 'g');
%             end
%             plot3(x_blade_plot_old, y_blade_plot_old, z_blade_plot_old, 'k');
%             hold off
%
%             view([0 90])
%
%             saveas(h,sprintf('./plot/frame/frame%d.png',i));


%             wake_panel_coord_old = wake_panel_coord;
%             [x_wake_panel_plot_old, y_wake_panel_plot_old, z_wake_panel_plot_old] = data2plot(wake_panel_coord_old(:, 1), wake_panel_coord_old(:, 2), wake_panel_coord_old(:, 3), 4);
%          blade_coord_old = blade_coord;
%         [x_blade_plot_old, y_blade_plot_old, z_blade_plot_old] = data2plot(blade_coord_old(:, 1), blade_coord_old(:, 2), blade_coord_old(:, 3), 4);

    end


        h = figure(1);
%         figure(1)
%         frame = plot3(x_blade_plot, y_blade_plot, z_blade_plot, 'b');
        plot3(x_blade_plot, y_blade_plot, z_blade_plot, 'g');
        if(~isempty(x_wake_panel_plot) && isempty(particles_coord))
            hold on
%             frame = plot3(x_wake_panel_plot, y_wake_panel_plot, z_wake_panel_plot, 'r');
            plot3(x_wake_panel_plot, y_wake_panel_plot, z_wake_panel_plot, 'k');
            hold off

        elseif(~isempty(particles_coord))
            hold on
%             frame = plot3(x_wake_panel_plot, y_wake_panel_plot, z_wake_panel_plot, 'r');
            plot3(x_wake_panel_plot, y_wake_panel_plot, z_wake_panel_plot, 'k');
            vec2plot = [];
            for j = 1:n_blades
                vec2plot = (((i - n_panels - 2) * j * (n_elements+1) - (i - n_panels - 3)):((i - n_panels - 2) * j * (n_elements+1) + (i - n_panels - 3) - (i - n_panels - 3)))';
                index = [];
%                 for k = 1:length(vec2plot)
%                     buono = find(abs(particles_coord(vec2plot(k), 3) - height_discretization) <= toll, 1);
%                     index(k) = buono;
%                 end
                scatter3(particles_coord(vec2plot, 1), particles_coord(vec2plot, 2), particles_coord(vec2plot, 3), 4, 'b', 'MarkerFaceColor',[0 .75 .75]);
            end
%             frame = scatter3(particles_coord(:, 1), particles_coord(:, 2), particles_coord(:, 3), 4);
%             scatter3(particles_coord(:, 1), particles_coord(:, 2), particles_coord(:, 3), 4, 'MarkerEdgeColor','k',...
%         'MarkerFaceColor',[0 .75 .75]);
%                coord1 = particles_coord(1:n_blades * (n_elements+1) * (i-3), :);
%                coord2 = particles_coord(n_blades * (n_elements+1) * (i-3) + 1:2*n_blades * (n_elements+1) * (i-3), :);
%                scatter3(coord1(:, 1), coord1(:, 2), coord1(:, 3), 4, 'MarkerEdgeColor','b',...
%           'MarkerFaceColor',[0 .75 .75]);
%                scatter3(coord2(:, 1), coord2(:, 2), coord2(:, 3), 4, 'MarkerEdgeColor','r',...
%           'MarkerFaceColor',[0 .75 .75]);
            hold off

        end

        axis([min_x max_x   min_y max_y   min_z max_z+0.05]);
        text(max_x-0.5, 1, max_z , t_plot_1)
        xlabel('x')
        ylabel('y')
        zlabel('z')
%          view([90 0])
%         view([90 -90])
         view([0 0])
%         view([-0.8355, -0.7, 1])
         saveas(h,sprintf('./plot/frame/frame%d.png',i));
% %
%         currFrame = getframe(gcf);
%
%         writeVideo(video, currFrame);
%
%         delete(frame)


end

% close(video);




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


%%
clc
figure(8)
for i = 1:n_timestep
    i
    blade_coord = importdata(strcat('./solution/blade_coord/coord', num2str(i), '.dat'));
    x_blade_plot = [];
    y_blade_plot = [];
    z_blade_plot = [];
    for j = 1:n_rotors
        coord = blade_coord((j-1) * n_blades * n_elements * 4 + 1: j *n_blades * n_elements * 4 , :);
        [x_blade_plot(:, j), y_blade_plot(:, j), z_blade_plot(:, j)] = data2plot(coord(:, 1), coord(:, 2), coord(:, 3), 4);
    end
    dFz = importdata(strcat('./solution/blade_Fz/Fz', num2str(i), '.dat'));
    blade_coord_Fz = importdata(strcat('./solution/blade_coord_Fz/coord', num2str(i), '.dat'));
    x_plot_Fz = reshape(blade_coord_Fz(:, 1), [n_elements, n_blades]);
    y_plot_Fz = reshape(blade_coord_Fz(:, 2), [n_elements, n_blades]);
    x_plot_Fz(end+1, :) = NaN;
    y_plot_Fz(end+1, :) = NaN;
    x_plot_Fz = x_plot_Fz(:);
    y_plot_Fz = y_plot_Fz(:);
    dFz(:, end + 1) = NaN;
    dFz = dFz';
    dFz2plot = dFz(:);
    plot3(x_plot_Fz, y_plot_Fz, dFz2plot)
    hold on
    plot3(x_blade_plot, y_blade_plot, z_blade_plot, 'g')
    axis([-1 1 -1 1 0 25])
    hold off
    pause(0.1)
end
