function cdf = empcdf_con(Data)
% Xu Shiyan, NTU
% convert pdf to its corresponding cdf


[sortedData oldIndex] = sort(Data);%sort in ascending order
cdf = ones(size(Data));
for i=1:size(cdf,2);%1~25
    m = 1;
    n = 1;
    for j=2:size(cdf,1)%2~250
          if sortedData(j,i)==sortedData(j-1,i)
              n = n + 1;
              continue;
          else
              for k=m:n
                  cdf(oldIndex(k,i),i) = j-1;
              end
              n = n + 1;
              m = n;
          end
    end
    for j=m:n
        cdf(oldIndex(j,i),i) = size(Data,1);
    end
end
cdf = (cdf-0.5)/size(sortedData,1);