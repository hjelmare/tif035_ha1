clear all
clc

a = [0.297104, 1.236745, 5.749982, 38.216677];

h = zeros(4,4);
for i = 1:4
    for j = 1:4
        prefactor = (2*pi)/(a(i) + a(j));
        pifactor = sqrt( pi/(a(i) + a(j)) );
        h(i,j) = prefactor * ( 3*a(j)*pifactor/2 - ...
                                3*a(j)^2*pifactor/(2*(a(i) + a(j)))  - ...
                                2 );
    end
end

Q = zeros(4,4,4,4);
for p = 1:4
    for r = 1:4
        for q = 1:4
            for s = 1:4
                Q(p,r,q,s) = 2*pi^(5/2) / ...
                    ( (a(p) + a(q))*(a(r) + a(s))*sqrt(a(p) + a(q) + a(r) + a(s)) );
            end
        end
    end
end


S = zeros(4,4);
for p = 1:4
    for q = 1:4
        S(p,q) = 4*pi*sqrt(pi/(a(p)+a(q)))/(4*(a(p)+a(q)));
    end
end

C = ones(4,1);
C = C/sqrt(C'*S*C);

F = zeros(4,4);

energyChange = 1;
oldEnergy = 10;
energy = 0;
while energyChange > 1e-8
    for p = 1:4
        for q = 1:4
            temp = 0;
            for r = 1:4
                for s = 1:4
                    temp = temp + Q(p,r,q,s)*C(r)*C(s);
                end
            end
            F(p,q) = h(p,q) + temp;
        end
    end

    [eigVectors,eigValues] = eig(F, S);

    eigValues = eigValues(eigValues ~= 0);
    [~,index] = min(eigValues);

    C = eigVectors(:,index);
    C = C/sqrt(C'*S*C);

    temp = 0;
    for p = 1:4
        for q = 1:4
            for r = 1:4
                for s = 1:4
                    temp = temp + Q(p,r,q,s)*C(p)*C(q)*C(r)*C(s);
                end
            end
        end
    end
    oldEnergy = energy;
    energy = 2*C'*h*C + temp;
    
    energyChange = abs(oldEnergy - energy);

end

energy

% Plotting stuff

nPoints = 1000;
rMax = 2;
radius = linspace(0,rMax,nPoints);

chi = @(rad) [exp(-a(1)*rad^2), exp(-a(2)*rad^2), exp(-a(3)*rad^2), exp(-a(4)*rad^2) ];
phi = @(rad) chi(rad)*C;
y = zeros(1,nPoints);
for ri = 1:nPoints
    y(ri) = phi(radius(ri))^2;
end

plot(radius, y)
xlabel('Radius [Å]','FontSize',14)
ylabel('\psi(r)','FontSize',14)
%saveas(gcf,'task1.png','png')


save('helium_coefficient.mat','C')
