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

tolerance = 1e-3;
oldEnergy = 1;
properEnergy = 2;

while abs(oldEnergy - properEnergy) > tolerance
    % Find Hartree potential (and Vxc) (task 2 + extension)
    
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
    % FUCKING MAGIC, HOW DOES THIS WORK? ;)
    U(1,1) = 1;
    U(1,2) = 0;
    U(end:end) = 1;
    U(end:end-1) = 0;
    u = U\(4*pi*radius.*n.*stepWidth^2)';

    % Translating back to reality
    V_sH = zeros(1,nPoints);
    for ri = 1:nPoints
        %V_sH(ri) = -2*u(ri)/radius(ri) + 1/rMax;   % THIS IS WHERE EXCHANGE HAPPENS??
        V_sH(ri) = -u(ri)/radius(ri) + 1/rMax;      % not sure about 2*
    end
    
    %---------------------------Enda skillnaden från task4---------------
    %Calculates V_c (eq 22-24)
    V_x = zeros(1,nPoints);
    e_x = zeros(1,nPoints);
    for ri = 1:nPoints
        e_x(ri) = -3/4*(3*n(ri)/pi)^(1/3);
        de_x = -(3/pi)^(1/3)/(4*n(ri)^(2/3));
        V_x(ri) = e_x(ri) + n(ri)*de_x;
    end
    %--------------------------------------------------------------------

    %--------------------new from task5-------------------
    [V_c, e_c] = GetV_c(n);
    %-----------------------------------------------------
    
    % Solve Kohn-Sham (task 3 + extension)
    H = zeros(nPoints);
    % Diagonal
    for i = 1:nPoints
        H(i,i) = 1/stepWidth^2 - 2/radius(i) + 2*V_sH(i) + V_x(i) + V_c(i);    % POTS GO HERE
    end

    % Sub- and superdiagonals
    for i = 2:nPoints
        H(i,i-1) = -.5/stepWidth^2;
        H(i-1,i) = -.5/stepWidth^2;
    end
    % TILLFÄLLIGT TILLAGDA MINUS ETT OCH SÅ
    H(1,1) = 1;
    H(1,2) = -0.00001;
    H(end,end)=1;
    H(end,end-1)=-0.00001;
    
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
    properEnergy = 2*gsEig - 2*trapz((.5*V_sH +V_x + V_c - e_x - e_c)'.*u.^2)*stepWidth;
    disp([norm, gsEig, peak, properEnergy])
    
    % Normalization - NOT SURE THIS IS VERY CLEVERLY DONE RIGHT NOW
    prefactor = 1/norm;
    n = prefactor*n;
        
    plot(radius,abs(gsWave)./radius')
    drawnow
end

disp('Final energy')
disp(properEnergy)




