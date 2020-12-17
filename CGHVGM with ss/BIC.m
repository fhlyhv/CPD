function [Kbic,nhbic] = BIC(GDat,l_array,g_array)

nl=size(l_array,2);
ng=size(g_array,2);
[n,p]=size(GDat);
K_array=zeros(p,p,ng*nl);
nh_array = zeros(ng*nl,1);
bic = zeros(ng*nl,1);
S=GDat'*GDat/n;

for j=1:nl
    for i=1:ng
        [X_1,X_2]=latentvariable(GDat',l_array(j),g_array(i),p);
        K=X_1+X_2;
        K(abs(K)<1e-4)=0;
        K_array(:,:,(j-1)*ng+i)=K;
        nh = rank(X_2);
        nh_array((j-1)*ng+i) = nh;
        bic((j-1)*ng+i) = -n*(logdet(X_1)-sum(sum(X_1.*S))) + log(n)*(nh*p+sum(sum(abs(sign(triu(K))))));
    end
end

Kbic = K_array(:,:,bic == min(bic));
nhbic = nh_array(bic==min(bic));