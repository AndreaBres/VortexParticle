close all
clear all
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
dt = 5e-4;
n_timestep = t_f / dt;
% t_f = 0.285;


%% Acquiring data 
particles_coord = [];

for i = 2:2:n_timestep
    i
        
        blade_coord = importdata(strcat('./solution/blade_coord/coord', num2str(i), '.dat'));
        [x_blade_plot, y_blade_plot, z_blade_plot] = data2plot(blade_coord(:, 1), blade_coord(:, 2), blade_coord(:, 3), 4);
        
        min_x = min(min(x_blade_plot));
        min_y = min(min(y_blade_plot));
        min_z = min(min(z_blade_plot));
        max_x = max(max(x_blade_plot));
        max_y = max(max(y_blade_plot));
        max_z = max(max(z_blade_plot));
        
    if i >= n_panels + 3 
       
        particles_coord = importdata(strcat('./solution/particle_intersect/coord', num2str(i), '.dat'));        
        
    end
    
    if(~isempty(particles_coord))
        min_x = min([min(x_blade_plot), min(particles_coord(:, 1)), min_x]);
        min_y = min([min(y_blade_plot), min(particles_coord(:, 2)), min_y]);
        min_z = min([min(z_blade_plot), min(particles_coord(:, 3)), min_z]);
        max_x = max([max(x_blade_plot), max(particles_coord(:, 1)), max_x]);
        max_y = max([max(y_blade_plot), max(particles_coord(:, 2)), max_y]);
        max_z = max([max(z_blade_plot), max(particles_coord(:, 3)), max_z]);
    end
         
end


%% Plotting hub

close all

width = 0.1;
offset_hub = 0.05;
R = 0.8;
R_circ = 0.99 * 0.25 * R;
dx = 0.001;
x = R_circ : -dx : -R_circ;
x(end+1:2*end-1) = -R_circ + dx : dx : R_circ;
circumference = @(x) sqrt(R_circ^2 - x.^2);

circ_up_coord = [];

circ_up_coord(1:length(x)/2) = circumference(x(1:end/2));
circ_up_coord(length(x)/2+1:length(x)) = - circumference(x(end/2+1:end));
x(end+1) = x(1);
circ_up_coord(end+1) = circ_up_coord(1);

z_circ_up = zeros(1, length(x)) - offset_hub;

x_vertexes =[R_circ, 0, -R_circ, 0];
y_vertexes =[0, -R_circ, 0, R_circ]; 
z_vertexes =[0, 0, 0, 0]- offset_hub;
% 
% figure
% patch(x, circ_up_coord, z_circ_up + width, 'red')
% hold on
% patch(x, circ_up_coord, z_circ_up - width, 'red')
% patch([x_vertexes(1), x_vertexes(2), x_vertexes(2), x_vertexes(1)],...
%       [y_vertexes(1), y_vertexes(2), y_vertexes(2), y_vertexes(1)],... 
%       [z_vertexes(1)+width, z_vertexes(2)+width, z_vertexes(2)-width, z_vertexes(1)-width], 'red')
% patch([x_vertexes(2), x_vertexes(3), x_vertexes(3), x_vertexes(2)],...
%       [y_vertexes(2), y_vertexes(3), y_vertexes(3), y_vertexes(2)],... 
%       [z_vertexes(2)+width, z_vertexes(3)+width, z_vertexes(3)-width, z_vertexes(3)-width], 'red')
% patch([x_vertexes(3), x_vertexes(4), x_vertexes(4), x_vertexes(3)],...
%       [y_vertexes(3), y_vertexes(4), y_vertexes(4), y_vertexes(3)],... 
%       [z_vertexes(3)+width, z_vertexes(4)+width, z_vertexes(4)-width, z_vertexes(3)-width], 'red')
% patch([x_vertexes(4), x_vertexes(1), x_vertexes(1), x_vertexes(4)],...
%       [y_vertexes(4), y_vertexes(1), y_vertexes(1), y_vertexes(4)],... 
%       [z_vertexes(4)+width, z_vertexes(1)+width, z_vertexes(1)-width, z_vertexes(4)-width], 'red')
% hold off

%% Plotting

particles_coord = [];
t_vec = dt:dt:t_f;

% video = VideoWriter('./plot/video/visual15.avi');
% video.FrameRate=15;
% open(video);

% for s = 1:2
s = 2;
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
    
    if i >= n_panels + 2
       
        particles_coord = importdata(strcat('./solution/particle_intersect/coord', num2str(i), '.dat'));     
    end
   
               
        h = figure(1);
%         figure(1)
%         frame = plot3(x_blade_plot, y_blade_plot, z_blade_plot, 'b');
        plot3(x_blade_plot, y_blade_plot, z_blade_plot, 'g');
            
        if(~isempty(particles_coord))
            hold on
           coord1 = particles_coord;
%            coord2 = particles_coord(n_blades * (n_elements+1) * (i-3) + 1:2*n_blades * (n_elements+1) * (i-3), :);
           scatter3(coord1(:, 1), coord1(:, 2), coord1(:, 3), 4, 'MarkerEdgeColor','b',...
          'MarkerFaceColor',[0 .75 .75]);
%                scatter3(coord2(:, 1), coord2(:, 2), coord2(:, 3), 4, 'MarkerEdgeColor','r',...
%           'MarkerFaceColor',[0 .75 .75]);
            hold off
            
        end
hold on
if s == 1
    % Hub plot 
    patch(x, circ_up_coord, z_circ_up + width, 'red')
    
    patch(x, circ_up_coord, z_circ_up - width, 'red')
end
    patch([x_vertexes(1), x_vertexes(2), x_vertexes(2), x_vertexes(1)],...
          [y_vertexes(1), y_vertexes(2), y_vertexes(2), y_vertexes(1)],... 
          [z_vertexes(1)+width, z_vertexes(2)+width, z_vertexes(2)-width, z_vertexes(1)-width], 'red')
    patch([x_vertexes(2), x_vertexes(3), x_vertexes(3), x_vertexes(2)],...
          [y_vertexes(2), y_vertexes(3), y_vertexes(3), y_vertexes(2)],... 
          [z_vertexes(2)+width, z_vertexes(3)+width, z_vertexes(3)-width, z_vertexes(3)-width], 'red')
    patch([x_vertexes(3), x_vertexes(4), x_vertexes(4), x_vertexes(3)],...
          [y_vertexes(3), y_vertexes(4), y_vertexes(4), y_vertexes(3)],... 
          [z_vertexes(3)+width, z_vertexes(4)+width, z_vertexes(4)-width, z_vertexes(3)-width], 'red')
    patch([x_vertexes(4), x_vertexes(1), x_vertexes(1), x_vertexes(4)],...
          [y_vertexes(4), y_vertexes(1), y_vertexes(1), y_vertexes(4)],... 
          [z_vertexes(4)+width, z_vertexes(1)+width, z_vertexes(1)-width, z_vertexes(4)-width], 'red')
    hold off


        
        axis([-1 1   -1 1   -1 1]);
% axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')

%          view([90 0])
%          view([90 -90])
if s == 1
%         view([-0.8355, -0.7, 1])
        text(1-0.2, 1-0.1, max_z+0.05, t_plot_1)
else
         view([0 90]) 
        text(1-0.2, 1-0.1, max_z+0.05, t_plot_1)
end
         saveas(h,sprintf('./plot/frame/frame%d.png',i)); 
% %         
%         currFrame = getframe(gcf);
%     
%         writeVideo(video, currFrame);
%     
%         delete(frame)


end

% close(video);



if s == 1
video = VideoWriter('./plot/video/visual_intersect3D.avi');
else 
video = VideoWriter('./plot/video/visual_intersect.avi');
end
video.FrameRate=15;
open(video);
for k = 2:2 : n_timestep
    
    if(exist(strcat('./plot/frame/frame', num2str(k), '.png')))
    
      filename = strcat('./plot/frame/frame', num2str(k), '.png');
      thisimage = imread(filename);
      writeVideo(video, thisimage);
      delete(filename)
    end
    
end

close(video);
% end


