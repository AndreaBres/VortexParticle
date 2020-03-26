close all; clear; clc

data = importdata('solution/flapping.dat');
dt = 5e-4;
time = dt:dt:0.3; % HARD CODED


%%
figure('name', 'Flapping angle')
subplot(4, 1, 1)
plot(time, rad2deg(data(1:3:end, 1)))
ylabel('$\beta_1 \  [deg]$', 'interpreter', 'latex')

subplot(4, 1, 2)
plot(time, rad2deg(data(1:3:end, 2)))
ylabel('$\beta_2 \  [deg]$', 'interpreter', 'latex')

subplot(4, 1, 3)
plot(time, rad2deg(data(1:3:end, 3)))
ylabel('$\beta_3 \  [deg]$', 'interpreter', 'latex')

subplot(4, 1, 4)
plot(time, rad2deg(data(1:3:end, 4)))
ylabel('$\beta_4 \  [deg]$', 'interpreter', 'latex')
xlabel('$time \ [s]$', 'interpreter', 'latex')



figure('name', 'Flapping angle all-in-one')
plot(time, rad2deg(data(1:3:end, 1)))
hold on;
plot(time, rad2deg(data(1:3:end, 2)))
plot(time, rad2deg(data(1:3:end, 3)))
plot(time, rad2deg(data(1:3:end, 4)))
ylabel('$\beta_4 \  [deg]$', 'interpreter', 'latex')
xlabel('$time \ [s]$', 'interpreter', 'latex')
hold off;

figure('name', 'Flapping velocity')
subplot(4, 1, 1)
plot(time, data(2:3:end, 1))
ylabel('$\dot{\beta}_1 \  [deg]$', 'interpreter', 'latex')

subplot(4, 1, 2)
plot(time, data(2:3:end, 2))
ylabel('$\dot{\beta}_2 \  [deg]$', 'interpreter', 'latex')

subplot(4, 1, 3)
plot(time, data(2:3:end, 3))
ylabel('$\dot{\beta}_3 \  [deg]$', 'interpreter', 'latex')

subplot(4, 1, 4)
plot(time, data(2:3:end, 4))
ylabel('$\dot{\beta}_4 \  [deg]$', 'interpreter', 'latex')
xlabel('$time \ [s]$', 'interpreter', 'latex')
