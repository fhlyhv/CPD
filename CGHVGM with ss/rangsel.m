function [l_r,g_r,K_array]=rangsel(M,l_lb,l_stepsize,l_ub,g_lb,g_stepsize,g_ub)

l=l_lb:l_stepsize:l_ub;
g=g_lb:g_stepsize:g_ub;
nl=size(l,2);
ng=size(g,2);
n=size(M,2);
z=zeros(ng,nl);
K_array=zeros(n,n,ng*nl);


for j=1:nl
    lambda=l(j);
    for i=1:ng
        gamma=g(i);
        [X_1,X_2]=latentvariable(M',lambda,gamma,n);
        S=X_1+X_2;
        S=cov_normalize(S^-1)^-1;
        S(abs(S)<=1e-4)=0;
        K_array(:,:,(j-1)*ng+i)=S;
        z(i,j)=sum(sum(abs(sign(S))))-n;
    end
end



i1=1;
j1=1;
while z(i1,j1)>=0.5*n*(n-1)
    i1=i1+1;
    j1=j1+1;
end
g_l=g(i1);
l_l=l(j1);

i2=ng;
j2=nl;
while z(i2,j2)==0;
    i2=i2-1;
    j2=j2-1;
end
g_u=g(i2);
l_u=l(j2);



l_r=[l_l,l_u];
g_r=[g_l,g_u];

