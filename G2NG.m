function NG_arti = G2NG(G_orig)

% Yu Hang, Xu Shiyan

cdfTable = normcdf(G_orig);

NG_arti = zeros(size(G_orig));
p = size(G_orig,2);
I=randperm(p);
for j=1:p
    switch mod(I(j),6)
		case 0
			NG_arti(:,j)=betainv(cdfTable(:,j),0.5,log(j));
		case 1
			NG_arti(:,j)=expinv(cdfTable(:,j),10);
		case 2
			NG_arti(:,j)=expinv(cdfTable(:,j),j);
		case 3
			NG_arti(:,j)=betainv(cdfTable(:,j),10*j,1);
		case 4
			NG_arti(:,j)=expinv(cdfTable(:,j),10*j);
		case 5
			NG_arti(:,j)=expinv(cdfTable(:,j),10+j,1);
    end
    if mod(j,5)==1
        NG_arti(:,j)=chi2inv(cdfTable(:,j),1);
    end
end