function [ V_c ] = GetV_c( density )
%Calculates the correlation part of the potential. requires the density (n)
%and returns the correlation potential. These calculations are based on 
% eq(22 - 27) in the task paper. 

gamma = -0.1423;
beta1 = 1.0529;
beta2 = 0.3334;
[~, width] = size(density);

%Calculates constants of the large expression that is de_c/dn
dt1 = beta1/(2*2^(1/3)*3^(5/6)*pi^(1/6));
dt2 = beta2/(6^(2/3)*pi^(1/3));
dn1 = (3/pi)^(1/6)*beta1/(2^(1/3));
dn2 = (3/pi)^(1/3)*beta2/(2^(2/3));

%constants for e_c
n1 = beta1*(3/(4*pi)^(1/6));
n2 = beta2*(3/(4*pi))^(1/3);

V_c = zeros(1, width);
for i = 1:width
    e_c = gamma/(1 + n1*(1/density(i))^(1/6) + n2*(1/density(i))^(1/3));
    
    de_c = gamma*(dt1*(1/density(i))^(7/6) + dt2*(1/density(i))^(4/3) )...
        /(dn1*(1/density(i))^(1/6) + dn2*(1/density(i))^(1/3) + 1);
    
    V_c(i) = e_c + density(i)*de_c;
    
    %Felsök, ta bort detta sen
    %---------------------------------------------------------
    r_s = (3/(4*pi*density(i)))^(1/3);
    if r_s < 1
        disp('Den blev mindre än 1------------------------------------')
    end
    %---------------------------------------------------------
end

end

