clear all
clc
clf

nPoints = 1000;
rMax = 10;

nIterations = 500000;

stepWidth = rMax / (nPoints-1);

radius = linspace(0,rMax,nPoints);

% Electron density - hydrogen
n = zeros(1,nPoints);
for ri = 1:nPoints
    n(ri) = exp(-2*radius(ri))/pi;
end

disp('Checking normalization - should be 1')
4*pi*trapz(n.*radius.*radius)*stepWidth

% Relaxation
U = zeros(nPoints);
u = zeros(nPoints,1);

for i = 1:nPoints
    U(i,i) = -2;
end
for i = 2:nPoints
    U(i-1,i) = 1;
    U(i,i-1) = 1;
end
% Boundary conditions
U(1,1) = 1;
U(1,2) = 0;
U(end:end) = 1;
U(end:end-1) = 0;
u = U\(4*pi*radius.*n.*stepWidth^2)';

% Translating back to reality
V = zeros(1,nPoints);
for ri = 1:nPoints
    V(ri) = -u(ri)/radius(ri) + 1/rMax;
end

% For reference / testing
V_hartree = zeros(1,nPoints);
for ri = 1:nPoints
    V_hartree(ri) = (1/radius(ri)) - (1 + (1/radius(ri)))*exp(-2*radius(ri));
end

% Plotting
hold on
plot(radius(1:10:end),V(1:10:end),'x')
plot(radius,V_hartree, 'r')
hold off
legend('Calculated potential using uniform grid', 'Analytically calculated', 'Density')
xlabel('Radius [Å]','FontSize',14)
ylabel('Hartree potential [au]', 'FontSize', 14)
%saveas(gcf,'task2_hydrogen.png','png')
