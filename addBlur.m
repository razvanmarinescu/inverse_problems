function imgNew = addBlur(img, blurSize)

% need to add blur too
imgNew = zeros(size(img));
[nrRows, nrCols] = size(img);

for r=1:nrRows
  for c=1:nrCols
    tmpSum = 0;
    count = 0;
    indRs = max(1,r-blurSize):min(nrRows,r+blurSize);
    indCs = max(1,c-blurSize):min(nrCols,c+blurSize);
    imgNew(r,c) = round(mean(mean(img(indRs, indCs))));
  end
end

end