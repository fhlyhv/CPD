function [X,pdf] = emppdf(XDat)

[cdf,X]=ecdf(XDat);

X=X(2:end);
pdf=cdf(2:end)-cdf(1:end-1);