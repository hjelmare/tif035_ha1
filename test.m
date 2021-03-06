clear all
clc

% Make initial guess at wave function (task 1)
task1



% and clean up stuff we won't need
clear C F Q S eigValues eigVectors energyChange h i index j oldEnergy
clear p q r s temp pifactor prefactor nPoints rMax radius ri y

hartreeToEV = 27.21;

nPoints = 1000;
rMin = 1e-10;
rMax = 5;
a0 = 1;
stepWidth = (rMax-rMin) / (nPoints-1);

radius = linspace(rMin,rMax,nPoints);

n = zeros(1,nPoints);
for ri = 1:nPoints
    n(ri) = phi(radius(ri))^2;
end

tolerance = 1e-4;
oldEnergy = 1;
properEnergy = 2;

u = sqrt(4*pi*n).*radius;
u = u/sqrt(trapz(radius, u.^2));

while abs(oldEnergy - properEnergy) > tolerance
    %V_sH = solveVSH(radius, u);
    V_sH = GetV_sH(u,radius);
    

    % Solve Kohn-Sham (task 3 + extension)
    H = zeros(nPoints);
    % Diagonal
    for i = 1:nPoints
        H(i,i) = 1/stepWidth^2 - 2/radius(i) + V_sH(i);    % POTS GO HERE
    end

    % Sub- and superdiagonals
    for i = 2:nPoints
        H(i,i-1) = -.5/stepWidth^2;
        H(i-1,i) = -.5/stepWidth^2;
    end
    % MINUS ETTOR WAAAT
    H(1,1) = 1;
    H(1,2) = 0;
    H(end,end)=1;
    H(end,end-1)=0;
    
    [vectors, values] = eig(H);

    
    
    values = sum(values);
    gsEig = min(values);
    gsIndex = find(values == gsEig);
    gsWave = vectors(:,gsIndex);
    
    for ri = 1:nPoints
        n(ri) = gsWave(ri)^2;
    end
    
    oldEnergy = properEnergy;
    % Checking/debugging inf
    disp('   Norm check         | eigenvalue         | peak of wf       | proper energy')
    norm = 4*pi*trapz(n.*radius.*radius)*stepWidth;
    peak = max(abs(gsWave));
    properEnergy = 2*gsEig - trapz(V_sH'.*u'.^2)*stepWidth;
    disp([norm, gsEig, peak, properEnergy])
    
    % Normalization
    prefactor = 1/norm;
    n = prefactor*n;
    
    
        
    plot(radius,abs(gsWave)./radius')
    drawnow
end
finalEnergy = properEnergy;
disp('Final energy')
disp(properEnergy)

