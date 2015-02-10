clear all
clc

rMin = 1e-3;
rMax = 1;
nPoints = 1000;
radius = linspace(rMin,rMax,nPoints);
stepWidth = rMax/(nPoints-1);

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

[vectors, values] = eig(H);

which = 2;

values = sort(unique(values));
values(which)

density = vectors(:,which);
for i = 1:nPoints
    density(i) = (density(i)/radius(i))^2 / (4*pi);
end

plot(radius,density)