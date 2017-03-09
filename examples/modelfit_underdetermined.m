x = [    87    91    93  ];%  92    93    95    98    99    99   105   105  101  80    89    91    92    93    94    95    99    99   104   105   110    86   90    91    92    93    95    97    99    99   105   105    94];
x = x';
ndat = size(x,1);
y = ones(ndat,1) + 2* x -1.3*sin(x.^2/30) + exp(-0.2*x.^3/40);
y = y.*(1 + 0.2*rand(ndat,1));
%y = [   76    71    40    43    69    84    38    41    30    22    38    37    90    42    39    43    51    79    33    33    50    30    32    24    59  31    70    67    59    76    30    27    23    41    35    53]';
[x2,ix] = sort(x);
y2 = y(ix);
xs = [min(x):max(x)];
figure;
hold on;
plot(x2,y2,'+');
for nfit = 2:2:6
    
    J = zeros(ndat,nfit);
    for k = 1:nfit
        J(:,k) = x2.^(k-1);
    end
    afit = J\y2;
    yfit = afit(1);
    for j = 2:nfit
        yfit = yfit + afit(j) * xs.^(j-1);
    end
    plot(xs,yfit,'r-');
end