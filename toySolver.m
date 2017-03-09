

%imgNewBlurred1d = reshape(addBlur(imgBlur, blurSize), [numel(imgBlur), 1 ]);
%solve the system (A^TA + alpha*I)f = A^Tg
%tol = 1e-6;
%maxit = 10;
lin1dFunc = @(x) [x(1) + 2*x(2), 2*x(1) + 4*x(2)]';

[fa,FLAG,RELRES,ITER,RESVEC] = pcg(lin1dFunc, [1,2]');


x1s = -2:0.1:2;
x2s = (1 - x1s)/2;

plot(x1s, x2s)
hold on
scatter(fa(1), fa(2))
xlim([0,2])
ylim([0,2])
