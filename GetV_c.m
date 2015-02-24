function [ V_c, e_c ] = GetV_c( density )
%Calculates the correlation part of the potential. requires the density (n)
%and returns the correlation potential. These calculations are based on 
% eq(22 - 27) in the task paper. 
A = 0.0311;
B = -0.048;
C = 0.0020;
D = -0.0116;
gamma = -0.1423;
beta1 = 1.0529;
beta2 = 0.3334;
[~, width] = size(density);

%constants for the case when r_s>=1-----------------------------
%Calculates constants of the large expression that is de_c/dn
dt1 = beta1/(2*2^(1/3)*3^(5/6)*pi^(1/6));
dt2 = beta2/(6^(2/3)*pi^(1/3));
dn1 = (3/pi)^(1/6)*beta1/(2^(1/3));
dn2 = (3/pi)^(1/3)*beta2/(2^(2/3));

%constants for e_c
n1 = beta1*(3/(4*pi)^(1/6));
n2 = beta2*(3/(4*pi))^(1/3);
%---------------------------------------------------------------

%Constants for the case when r_s <1------
%constants for de_c
dT1 = 6^(2/3)*pi^(1/3)*A;
dT2 = 3 - log(pi) + log(3) - 2*log(2);
dN1 = 3*6^(2/3)*pi^(1/3);

%constants for e_c =
constant1 = (3/(4*pi))^(1/3);
constant2 = C*constant1;
constant3 = D*constant1;
%----------------------------------------

V_c = zeros(1, width);
e_c = zeros(1, width);
for i = 1:width
    r_s = (3/(4*pi*density(i)))^(1/3);
    
    if r_s >=1
        e_c(i) = gamma/(1 + n1*(1/density(i))^(1/6) + n2*(1/density(i))^(1/3));
        
        de_c = gamma*(dt1*(1/density(i))^(7/6) + dt2*(1/density(i))^(4/3) )...
            /(dn1*(1/density(i))^(1/6) + dn2*(1/density(i))^(1/3) + 1);
        
        V_c(i) = e_c(i) + density(i)*de_c;
    else
        e_c(i) = A*log(constant1*(1/density(i))^(1/3)) + B + constant2*...
            log(constant1*(1/density(i))^(1/3)) + constant3*(1/density(i))^(1/3);
        
        de_c = -(dT1 + (1/density(i))^(1/3)*(C*(log(1/density(i)) + dT2) + 3*D)...
            )/(dN1*density(i));
        
        V_c(i) = e_c(i) + density(i)*de_c; 
    end
end
V_c(1) = 0;
V_c(end) = 0;

end

