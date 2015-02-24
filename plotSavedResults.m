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

load('task1_wave.mat')
phi1 = phiValues;
load('radius1.mat')
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



phi4 = u4./(sqrt(4*pi)*radius);
phi5 = u5./(sqrt(4*pi)*radius);
phi6 = u6./(sqrt(4*pi)*radius);

hold on
plot(radius1, phi1, 'g')
plot(radius, phi4, 'k:')
plot(radius, phi5, 'r')
plot(radius, phi6, 'b-.')

%plot(radius, u4, 'k:') %The radial wave functions (without task 1)
%plot(radius, u5, 'r')
%plot(radius, u6, 'b-.')

xlabel('radius [au]', 'FontSize', 14)
ylabel('u(r)', 'FontSize', 14)
legend('task1', 'task4', 'task5', 'task6')

