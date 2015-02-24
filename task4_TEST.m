clear all
clc

% Make initial guess at wave function (task 1)
task1

% and clean up stuff we won't need
clear C F Q S eigValues eigVectors energyChange h i index j oldEnergy
clear p q r s temp pifactor prefactor nPoints rMax radius ri y

% Convergence testing stuff
rMaxes = [3:0.2:10];

for outerIterations = 1:length(rMaxes)

    rMin = 1e-10;
    stepWidth = 0.05;
    rMax = rMaxes(outerIterations);
    nPoints = (rMax-rMin)/stepWidth + 1;
    nPoints = round(nPoints);
 
    radius = linspace(rMin,rMax,nPoints);
    % Get density from guesstimated wave function
    n = zeros(1,nPoints);
    for ri = 1:nPoints
        n(ri) = phi(radius(ri))^2;
    end
    % and then find u from the density
    u = sqrt(4*pi*n).*radius;
    u = u/sqrt(trapz(radius, u.^2));

    tolerance = 1e-5;   % Termination criterion
    oldEnergy = 1;      % Just some energies to get the loop started
    properEnergy = 2;

    while abs(oldEnergy - properEnergy) > tolerance
        V_sH = GetV_sH(u, radius);

        % Solve Kohn-Sham (task 3 + extension)
        H = zeros(nPoints);
        % Diagonal
        for i = 1:nPoints
            H(i,i) = 1/stepWidth^2 - 2/radius(i) + V_sH(i); % add pots here
        end

        % Sub- and superdiagonals
        b = -.5/stepWidth^2;
        for i = 2:nPoints
            H(i,i-1) = b;
            H(i-1,i) = b;
        end
        % Boundary conditions
        H(1,1) = 1;
        H(1,2) = 0;
        H(end,end)=1;
        H(end,end-1)=0;
        
        % Solve eigenproblem
        [vectors, values] = eig(H);
        
        %Get ground state
        values = diag(values);
        gsEig = min(values);
        gsIndex = find(values == gsEig);
        gsWave = vectors(:,gsIndex);

        u = gsWave';
        u = u/sqrt(trapz(radius, u.^2));

        oldEnergy = properEnergy;
        % Checking/debugging inf
        disp('   eigenvalue         | peak of wf       | proper energy')
        peak = max(abs(gsWave));
        properEnergy = 2*gsEig - trapz(V_sH.*u.^2)*stepWidth;
        disp([gsEig, peak, properEnergy])

    end
    finalEnergy = properEnergy;
    disp('Final energy')
    disp(properEnergy)
    
    %plotEnergy(outerIterations) = properEnergy;

end

%save('task4_EnergyOver_stepWidth.m', 'plotEnergy')
%save('stepWidth.m', 'stepWidths')
% % Plot for rMax sweep
%plot(rMaxes,plotEnergy)
%axis([5 15 -2.8535 -2.8531])
%xlabel('rMax [au]','FontSize',14)

% Plot for stepWidth sweep
%plot(stepWidths,plotEnergy)
%xlabel('stepWidth [au]','FontSize',14)


%ylabel('Energy [au]','FontSize',14)
