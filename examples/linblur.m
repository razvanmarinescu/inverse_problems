function [y,K] = linblur(f,s)
n = size(f,1);
K = zeros(n,n);
h = 1/n;
C = h/sqrt(2*pi*s*s);
x = [1:n];
k = C*exp(-x.^2 * h*h/(2*s*s));
for i = 1:n
    K(i,1:i) = k(i:-1:1);
    K(i,i+1:n) = k(2:n-i+1);
end
y = K*f;