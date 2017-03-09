function [Df] = lindf(f)
n = size(f,1);
Df = zeros(n+1,n);
h = 1/n;
x = [1:n];
k = [1,-1,zeros(1,n-2)];
for i = 1:n
    Df(i,1:i) = k(i:-1:1);
%    Df(i,i+1:n) = k(2:n-i+1);
end
Df(n+1,n) =1;