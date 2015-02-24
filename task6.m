clear all
clc

% Make initial guess at wave function (task 1)
task1

% and clean up stuff we won't need
clear C F Q S eigValues eigVectors energyChange h i index j oldEnergy
clear p q r s temp pifactor prefactor nPoints rMax radius ri y

rMin = 1e-10;
rMax = 8;
stepWidth = 0.004;
nPoints = (rMax-rMin)/stepWidth + 1;
nPoints = round(nPoints)

radius = linspace(rMin,rMax,nPoints);

n = zeros(1,nPoints);
for ri = 1:nPoints
    n(ri) = phi(radius(ri))^2;
end

u = sqrt(4*pi*n).*radius;
u = u/sqrt(trapz(radius, u.^2));

tolerance = 1e-5;
oldEnergy = 1;
properEnergy = 2;

while abs(oldEnergy - properEnergy) > tolerance
    V_sH = GetV_sH(u, radius);
    
    %Calculates V_c (eq 22-24)
    n = 2*((u./radius).^2/(4*pi));
    V_x = -(3.*n ./ pi).^(1/3);
    e_x = -3/4*(3*n/pi).^(1/3);
    
    % new from task5 
    [V_c, e_c] = GetV_c(n);
    
    % Solve Kohn-Sham (task 4 + extension)
    H = zeros(nPoints);
    % Diagonal
    for i = 1:nPoints
        H(i,i) = 1/stepWidth^2 - 2/radius(i) + 2*V_sH(i) + V_x(i) + V_c(i);
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
    
    u = gsWave';
    u = u/sqrt(trapz(radius, u.^2));
        
    oldEnergy = properEnergy;
    % Checking/debugging inf
    disp('   eigenvalue         | peak of wf       | proper energy')
    peak = max(abs(gsWave));
    properEnergy = 2*gsEig - 2*trapz((0.5*2*V_sH + V_x + V_c - e_x - e_c).*u.^2)*stepWidth;
    disp([gsEig, peak, properEnergy])
       
%     %plot(radius,abs(gsWave)./radius')
%     plot(radius,u.^2 ./ radius.^2 ,'r',radius,4*pi*n,'b')
%     drawnow
end
finalEnergy = properEnergy;
disp('Final energy')
disp(properEnergy)

