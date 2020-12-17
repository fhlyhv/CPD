% Xu Shiyan, 09/27/2011
% It reads non-Gaussian data and outputs standard Gaussian data

function stdNormData = copula(rawData)
[sortedData oldIndex] = sort(rawData);
cdfTable = zeros(size(rawData));
for i=1:size(cdfTable,2)
    m = 1;
    n = 1;
    for j=2:size(cdfTable,1)
          if sortedData(j,i)==sortedData(j-1,i)
              n = n + 1;
              continue;
          else
              for k=m:n
                  cdfTable(oldIndex(k,i),i) = j-1;
              end
              n = n + 1;
              m = n;
          end
    end
    for j=m:n
        cdfTable(oldIndex(j,i),i) = size(rawData,1)-1e-8;
    end
end
cdfTable = cdfTable./size(sortedData,1);
stdNormData = norminv(cdfTable,0,1);
end