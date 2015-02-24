clear all
clc

% Make initial guess at wave function (task 1)
task1



% and clean up stuff we won't need
clear C F Q S eigValues eigVectors energyChange h i index j oldEnergy
clear p q r s temp pifactor prefactor nPoints rMax radius ri y

hartreeToEV = 27.21;

stepWidth = 5e-3;
rMin = 1e-10;
rMax = 12;


nPoints = (rMax-rMin) / stepWidth + 1;
nPoints = round(nPoints)

radius = linspace(rMin,rMax,nPoints);

n = zeros(1,nPoints);
for ri = 1:nPoints
    n(ri) = phi(radius(ri))^2;
end

tolerance = 1e-6;
oldEnergy = 1;
properEnergy = 2;

while abs(oldEnergy - properEnergy) > tolerance
    % Find Hartree potential (and Vxc) (task 2 + extension)
    
    % Relaxation
    operatorU = zeros(nPoints);
    U = zeros(nPoints,1);
    
    for i = 1:nPoints
        operatorU(i,i) = -2;
    end
    for i = 2:nPoints
        operatorU(i-1,i) = 1;
        operatorU(i,i-1) = 1;
    end
    % FUCKING MAGIC, HOW DOES THIS WORK? ;)
    operatorU(1,1) = 1;
    operatorU(1,2) = 0;
    operatorU(end:end) = 1;
    operatorU(end:end-1) = 0;
    U = operatorU\(4*pi*radius.*n.*stepWidth^2)';
    
    % Translating back to reality
    V_sH = zeros(1,nPoints);
    for ri = 1:nPoints
        V_sH(ri) = -U(ri)/radius(ri) + 1/rMax;
    end

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
    
    H(1,1) = 1;
    H(1,2) = 0;
    H(end,end)=1;
    H(end,end-1)= 0;

    [vectors, values] = eig(H);
    
    values = diag(values);
    gsEig = min(values);
    gsIndex = find(values == gsEig);
    gsWave = vectors(:,gsIndex);
    
    for ri = 1:nPoints
        n(ri) = gsWave(ri)^2 / (4*pi*radius(ri)^2);
    end
        
    oldEnergy = properEnergy;
    
    % Checking/debugging inf
    disp('   Norm check         | eigenvalue         | peak of wf       | proper energy')
    norm = 4*pi*trapz(n.*radius.*radius)*stepWidth;
    prefactor = 1/norm;
    n = prefactor*n;
    peak = max(abs(gsWave));
    properEnergy = 2*gsEig - trapz(V_sH'.*(sqrt(4*pi*n').*radius').^2)*stepWidth;
    disp([norm, gsEig, peak, properEnergy])
        
    plot(radius,abs(gsWave)./radius')
    drawnow
end
properEnergy = 2*gsEig - trapz(V_sH'.*(sqrt(4*pi*n').*radius').^2)*stepWidth;

disp('Final energy')
disp(properEnergy)
