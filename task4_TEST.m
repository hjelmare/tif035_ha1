clear all
clc

% Make initial guess at wave function (task 1)
task1



% and clean up stuff we won't need
clear C F Q S eigValues eigVectors energyChange h i index j oldEnergy
clear p q r s temp pifactor prefactor nPoints rMax radius ri y


%rMin = 1e-10;
%rMax = 6;
%stepWidth = 0.001;
%nPoints = (rMax-rMin)/stepWidth + 1;
%nPoints = round(nPoints)

%rMaxes = [5:15];
stepWidths = [0.1 0.05 0.04 0.02 0.01 0.005 0.004 0.003 0.001 0.0005];



for outerIterations = 1:length(stepWidths)

    rMin = 1e-10;
    %stepWidth = 0.05;
    rMax = 6;
    %rMax = rMaxes(outerIterations);
    stepWidth = stepWidths(outerIterations);
    nPoints = (rMax-rMin)/stepWidth + 1;
    nPoints = round(nPoints);
    
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
        %V_sH = solveVSH(radius, u);

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
        %disp('   Norm check         | eigenvalue         | peak of wf       | proper energy')
        %norm = 4*pi*trapz(n.*radius.*radius)*stepWidth;
        peak = max(abs(gsWave));
        properEnergy = 2*gsEig - trapz(V_sH.*u.^2)*stepWidth;
        %disp([0000, gsEig, peak, properEnergy])

    end
    finalEnergy = properEnergy;
    %disp('Final energy')
    %disp(properEnergy)
    
    plotEnergy(outerIterations) = properEnergy;

end
save('task4_EnergyOver_stepWidth.m', 'plotEnergy')
save('stepWidth.m', 'stepWidths')
% % Plot for rMax sweep
%plot(rMaxes,plotEnergy)
%axis([5 15 -2.8535 -2.8531])
%xlabel('rMax [au]','FontSize',14)

% Plot for stepWidth sweep
plot(stepWidths,plotEnergy)
xlabel('stepWidth [au]','FontSize',14)


ylabel('Energy [au]','FontSize',14)
