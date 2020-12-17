function lambda = lambda_sel(K_est,K_array,l_lb,l_stepsize,l_ub,g_lb,g_stepsize,g_ub)

l=l_lb:l_stepsize:l_ub;
g=g_lb:g_stepsize:g_ub;
nl=size(l,2);
ng=size(g,2);
hd=zeros(1,nl*ng);


for i=1:nl*ng
    hd(i)=sum(sum(abs(abs(sign(K_array(:,:,i)))-K_est)));
end

hd=reshape(hd,ng,nl);

[~,b]=find(hd==min(hd(:)));

lambda=mean(l(b));

