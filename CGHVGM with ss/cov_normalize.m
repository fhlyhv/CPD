function SS=cov_normalize(S)

n=size(S,1);
A=ones(1,n);
for i=1:n
    A(i)=S(i,i);
end
SS=zeros(n);
for i=1:n
    for j=1:n
        SS(i,j)=S(i,j)/(sqrt(A(i)*A(j)));
    end
end
    