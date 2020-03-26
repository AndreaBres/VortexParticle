close all
clear all
clc

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
dt = 5e-4;
R = 0.8;
% t_f = 0.285;

n_timesteps = t_f / dt;

Re = importdata('./solution/Re.dat');

t_vec = dt:dt:t_f;

figure(1)
for i = 1:n_rotors
    Re_vec(i, :) = Re((i-1) * n_timesteps + 1 : i * n_timesteps);
end
for i = 1:n_rotors
    plot(t_vec, Re_vec(i, :)/R)
    hold on
end

title('$R_e$ nel tempo della simulazione', 'Interpreter', 'latex')
xlabel('time [s]', 'Interpreter', 'latex')
ylabel('$R_e$/R', 'Interpreter', 'latex')

