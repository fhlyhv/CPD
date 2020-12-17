function C = Contra(GDat,N)
%The contrast function
%Li Chenyang, NTU, Jun. 2013

%#codegen

[n,p]=size(GDat);
aa = GDat'*GDat/n;

while 1
    [~,q] = chol(aa);
    if q == 0
        break;
    else
        aa = aa+1e-4*speye(p);
    end
end


C =n/N*2*sum(log(diag(chol(aa))));  %logdet(GDat.'*GDat/n);    %log(det(GDat.'*GDat/n)); %


