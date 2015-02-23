%Detta baseras på ekvationerna 22-27 i uppgiftsbladet

A = 0.0311;
B = -0.048;
C = 0.0020;
D = -0.0116;
gamma = -0.1423;
beta1 = 1.0529;
beta2 = 0.3334;

%the variable r_s is derived from the density (n);
n = 0.2; %change this!
r_s = (3/(4*pi*n))^(1/3); %??????????????? stämmer detta

e_c =gamma/(1+beta1*(3/(4*pi*n))^(1/6) + beta2*(3/(4*pi*n))^(1/3));

%OBS kontrollera tecken
%de_c = - gamma*(- beta1*(1/n)^(7/6)/(2*2^(1/3)*3^(5/6)*pi^(1/6)) - beta2*(1/n)^(4/3)/(6^(2/3)*pi^(1/3)) ) ...
%        /((3/pi)^(1/6)*beta1*(1/n)^(1/6)/(2^(1/3)) + (3/pi)^(1/3)*beta2*(1/n)^(1/3)/(2^(2/3)) + 1)^2;
    
de_c = gamma*(beta1*(1/n)^(7/6)/(2*2^(1/3)*3^(5/6)*pi^(1/6)) + beta2*(1/n)^(4/3)/(6^(2/3)*pi^(1/3)) ) ...
        /((3/pi)^(1/6)*beta1*(1/n)^(1/6)/(2^(1/3)) + (3/pi)^(1/3)*beta2*(1/n)^(1/3)/(2^(2/3)) + 1)^2;
    
%om r_s < 1
e_c = A*ln((3/(4*pi*x))^(1/3)) + B + C*(3/(4*pi*x))^(1/3)*ln((3/(4*pi*x))^(1/3)) + D*(3/(4*pi*x))^(1/3);
de_c = DU HAR SKRIVIT NER VAD DET SKA VARA I DITT BLOCK!!!!!!!1
