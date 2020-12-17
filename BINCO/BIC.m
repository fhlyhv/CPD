function Kbic=BIC(XDat,l_array)

[n,p] = size(XDat);
S = XDat.'*XDat/n;
nl = length(l_array);

bic = zeros(nl,1);
Ktmp = zeros(p,p,nl);

for k = 1:nl
    K = L1precisionBCD(S,exp(l_array(k)));
    K(abs(K)<1e-4)=0;
    bic(k) = -n*(logdet(K)-sum(sum(sparse(K).*S)))+log(n)*sum(sum(abs(sign(triu(K))))); 
    Ktmp(:,:,k) = K;
end

Kbic = Ktmp(:,:,bic == min(bic));