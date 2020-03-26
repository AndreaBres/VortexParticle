

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
t_f = 0.25;
dt = 5e-4;
Vtip = 100.53;
R = 0.8;

n_timesteps = t_f / dt;

%%

max_vi = -100;
min_vi = 100;

for i = 1:n_timesteps
    
    filename = strcat('./solution/vi_on_blade/vi_on_blade', num2str(i), '.dat');
    vi_on_blade = importdata(filename);
    max_vi = max(max_vi, max(vi_on_blade(:, 3)));
    min_vi = min(min_vi, min(vi_on_blade(:, 3)));
    
end

%% Plotting induced velocity

Omega = rad2deg(Vtip / R);
time_vec = dt:dt:t_f;
psi = 0:360/n_blades:(360-360/n_blades);


for i = 1:n_timesteps
% i = 1;
    psi = psi + Omega * dt;
    psi(psi > 360) = psi(psi > 360) - 360;
    % Acquiring datas of vi 
    filename = strcat('./solution/vi_on_blade/vi_on_blade', num2str(i), '.dat');
    vi_on_blade = importdata(filename);
    % Reshaping into a more useful matrix
    vi_on_blade = reshape(vi_on_blade, [n_elements+1, n_blades, 3]);
    x_plot = linspace(0.2, 0.8, n_elements+1);
    h = figure(1);
    sgtitle(strcat('Induced velocity normal to the blades', {' '}, 't = ', num2str(time_vec(i))));
    for j = 1:n_blades
        % Required that n_blades is even
        k = double(j);
        for s = 1:n_elements + 1
            vi_2plot(s) = vi_on_blade(s, j, 3);
        end
        subplot(2, 2, k)
        plot(x_plot, vi_2plot)
        xlabel('x [m]')
        ylabel('vi [m/s]')
        axis([min(x_plot) max(x_plot) min_vi max_vi]);
        title(strcat('blade n°', num2str(j), {' '}, 'psi = ', num2str(psi(j)), '°'))
    end
       saveas(h,sprintf('./plot/frame/frame_vi/frame%d.png',i)); 
       close(h);
    
end

video = VideoWriter('./plot/video/visual_vi.avi');
video.FrameRate=15;
open(video);
for j = 1 : n_timesteps
    
    if(exist(strcat('./plot/frame/frame_vi/frame', num2str(j), '.png')))
    
      filename = strcat('./plot/frame/frame_vi/frame', num2str(j), '.png');
      thisimage = imread(filename);
      writeVideo(video, thisimage);
      delete(filename)
    end
    
end

close(video);
