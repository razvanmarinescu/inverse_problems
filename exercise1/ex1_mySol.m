

phi = @(x,p) sum(abs(x) .^ p);

A = [1, 2];
b = 5;

for p=1:0.5:4
  funcP = @(x) phi(x,p);

  [x, fVal] = fmincon(funcP, [0, 0]', [], [], A, b);

end
x2 = A\b;

%% Problem 2
n = 10;
xs = linspace(0, (2*n-1), 2*n)/(2*n-1);
ys = linspace(0, n-1, n)/(n-1);
sigma = 0.2;
matA = zeros(n, 2*n);
for i=1:n
  for j=1:(2*n)
    matA(i,j) = ys(i)-xs(j);
  end
end

matA = 1/(sqrt(2*pi)*sigma)*...
  exp(-(matA .^ 2)/(2*sigma^2));

imagesc(matA);

Aimg = ceil(matA/max(matA(:))*256);
colormap = parula(256);
imwrite(Aimg, colormap, 'Aimage2.png');

[v,w,u] = svd(matA);
[wSparse, d] = spdiags(pinv(w));

pinvA = pinv(matA)






