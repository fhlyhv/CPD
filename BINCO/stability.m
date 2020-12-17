function Kstb=stability(S,lambda,nid)


K_tmp= zeros(size(S));

for k=1:nid
    K=L1precisionBCD(S(:,:,k),lambda);
    K=inv(cov_normalize(inv(K)));
    K(abs(K)<1e-4)=0;
    K_tmp(:,:,k)=abs(sign(K));
end

Kstb=sum(K_tmp,3)/nid;