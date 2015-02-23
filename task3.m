clear all
clc

rMin = 1e-10;
rMax = 5;

nPoints = 1000;
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

H(1,1) = 1;
H(1,2) = 0;
H(end,end)=1;
H(end,end-1)=0;

[vectors, values] = eig(H);

values = sum(values);
gsEig = min(values);
gsIndex = find(values == gsEig);
gsWave = vectors(:,gsIndex);

plot(radius, abs(gsWave)./(sqrt(4*pi)*radius'))
disp('Ground state energy [eV]')
disp(gsEig*27.21)