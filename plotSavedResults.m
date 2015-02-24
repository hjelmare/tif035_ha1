%% Plot the energydependence of the step width (i.e. check convergense)
clear all
clc

load('dataTask4_energy.mat')
load('dataTask4_stepWidth.mat')
plot(b, a)

xlabel('Step width [au]', 'FontSize', 14)
ylabel('Energy [au]', 'FontSize', 14)
legend('Energy as a function of the step width')

%% Plot the resulting wave functions from task 4, 5 & 6
clear all
clc
close all

load('task4_wave.mat', 'u');
u4 = u;
clear 'u'
load('task5_wave.mat', 'u');
u5 = u;
clear 'u'
load('task6_wave.mat', 'u');
u6 = u;
clear 'u'
load('radius.mat')

hold on
plot(radius, u4, 'k:')
plot(radius, u5, 'r')
plot(radius, u6, 'b-.')

xlabel('radius [au]', 'FontSize', 14)
ylabel('u(r)', 'FontSize', 14)
legend('task4', 'task5', 'task6')

