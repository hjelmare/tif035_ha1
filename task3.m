clear all
clc
close all

rMin = 1e-10;
rMax = 8;

nPoints = 2000;
radius = linspace(rMin,rMax,nPoints);
stepWidth = (rMax-rMin)/(nPoints-1);

H = zeros(nPoints);
% Diagonal
for i = 1:nPoints
    H(i,i) = 1/stepWidth^2 - 1/radius(i);
end

% Sub- and superdiagonals
for i = 2:nPoints
    H(i,i-1) = -.5/stepWidth^2;
    H(i-1,i) = -.5/stepWidth^2;
end

% Boundary conditions
H(1,1) = 1;
H(1,2) = 0;
H(end,end)=1;
H(end,end-1)=0;

% Solve eigenproblem
[vectors, values] = eig(H);

values = sum(values);
gsEig = min(values);
gsIndex = find(values == gsEig);
gsWave = vectors(:,gsIndex);

HydrogenWaveFunc = @(r) r*exp(-r); % Analytic result from course book
u_H = zeros(1,length(radius));
for i = 1:length(radius)
    u_H(i) = HydrogenWaveFunc(radius(i));
end

%Normalization
u_H = u_H/sqrt(trapz(radius, u_H.^2));
gsWave = gsWave/sqrt(trapz(radius, gsWave.^2));

hold on
plot(radius, gsWave, 'o')
plot(radius, u_H, 'rd', 'MarkerSize', 2)
xlabel('radius [au]', 'FontSize', 14)
ylabel('u(r) [au]', 'FontSize', 14)
text = legend('Result from simulation', 'Analytically calculated');
set(text, 'FontSize', 12)

disp('Ground state energy [eV]')
disp(gsEig*27.21)


