function K_est=stasel_rnran(M,subn,alpha1,alpha2,lb_l,stepsize_l,ub_l,lb_g,stepsize_g,ub_g,EV)


[m,n]=size(M);
a=floor(m/2);
MM=zeros(a,n,subn);
KK=zeros(n,n,subn);
for kk=1:2:subn-1
    I=randperm(m);
    MM(:,:,kk)=M(I(1:a),:);
    MM(:,:,kk+1)=M(I(a+1:2*a),:);
end
lambda=lb_l:stepsize_l:ub_l;
gamma=lb_g:stepsize_g:ub_g;
nl=size(lambda,2);
ng=size(gamma,2);
K_ave=zeros(n,n,nl*ng);
q=zeros(1,nl*ng);




for j=1:nl
    l=lambda(j);
    for i=1:ng
        g=gamma(i);
        [X1q,X2q]=latentvariable(M',l,g,n);
        S=X1q+X2q;
        S=cov_normalize(S^-1)^-1;
        S(abs(S) < 1e-4) = 0;
        q((j-1)*ng+i)=(sum(sum(abs(sign(S))))-n)/2;
        
        rand('state',sum(100*clock)*rand(1));
        rl=rand(1,subn);
        rg=rand(1,subn);
        wl=alpha1+rl*(1-alpha1);
        wg=alpha2+rg*(1-alpha2);

        for kk=1:subn
            [X1,X2]=latentvariable(MM(:,:,kk)',l/wl(kk),g/wg(kk),n);
            K=X1+X2;
            K=cov_normalize(K^-1)^-1;
            K(abs(K) < 1e-4) = 0;
            KK(:,:,kk)=abs(sign(K));
        end
        K_ave(:,:,(j-1)*ng+i)=sum(KK,3)/subn;
 
    end
end

Eq=sum(q)/(ng*nl); 
pithr=((Eq^2)/(n*(n-1)*EV))+0.5;

K_ave(K_ave<=pithr)=0;
K_ave(K_ave>pithr)=1;
K_est=sign(sum(K_ave,3));

