function [ V_sH ] = GetV_sH( u, r)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
rMax = r(end);
[~, nPoints] = size(r);
h = r(3)-r(2);

U = zeros(nPoints);
for i = 1:nPoints
    U(i,i) = -2/(h^2);
end
for i = 2:nPoints
    U(i-1,i) = 1/(h^2);
    U(i,i-1) = 1/(h^2);
end

U(1,1) = 1;
U(1,2) = 0;
U(end,end) = 1;
U(end,end-1) = 0;

B = -u.^2./r;
A = U\B';

%Convert to V_sH
V_sH = A'./r + 1/rMax;
end

