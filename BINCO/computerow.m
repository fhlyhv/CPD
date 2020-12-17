function transCost = computerow(i,GDat,ndefault,cj,lt,n)
transCost = Inf*ones(1,lt);
ck = cj;
j = cj+ndefault*(i-1);
for k = j+ndefault-1:n-ndefault-1
    transCost(ck) = Contra(n,GDat(j:k,:));
    ck = ck+1;
end