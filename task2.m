clear all
clc
clf

nPoints = 1000;
rMax = 10;
a0 = 0.529;     % ångström??
a0 = 1;

nIterations = 500000;

stepWidth = rMax / (nPoints-1);

radius = linspace(0,rMax,nPoints);

% Electron density - hydrogen
n = zeros(1,nPoints);
for ri = 1:nPoints
    n(ri) = exp(-2*radius(ri))/pi;
end

% % Electron density - helium
% load('helium_coefficient.mat')
% a = [0.297104, 1.236745, 5.749982, 38.216677];
% chi = @(rad) [exp(-a(1)*rad^2), exp(-a(2)*rad^2), exp(-a(3)*rad^2), exp(-a(4)*rad^2) ];
% phi = @(rad) chi(rad)*C;
% n = zeros(1,nPoints);
% for ri = 1:nPoints
%     n(ri) = phi(radius(ri))^2;
% end

disp('Checking normalization - should be 1')
4*pi*trapz(n.*radius.*radius)*stepWidth

% Relaxation
U = zeros(1,nPoints);

U(1) = 0;
U(end) = 0;

for i = 1:nIterations
    for ri = 1+1 : nPoints - 1
        U(ri) = (1/2)*(4*pi*radius(ri)*n(ri)*stepWidth^2 + U(ri+1) + U(ri-1));
    end
end

% Translating back to reality
V = zeros(1,nPoints);
for ri = 1:nPoints
    V(ri) = U(ri)/radius(ri) + 1/rMax;
end

% For reference / testing
V_hartree = zeros(1,nPoints);
for ri = 1:nPoints
    V_hartree(ri) = (1/radius(ri)) - (1 + (1/radius(ri)))*exp(-2*radius(ri));
end


% Plotting
hold on
plot(radius,V)
plot(radius,V_hartree, 'r')
plot(radius,n,'k')
hold off
legend('Calculated potential', 'Hartree potential', 'Density')
xlabel('Radius [Å]','FontSize',14)
%saveas(gcf,'task2_hydrogen.png','png')
