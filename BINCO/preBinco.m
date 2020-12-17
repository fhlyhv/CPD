
tao = unique([1,tao,n]);
nt = length(tao)-1;
Kg = zeros(p,p,nt); % glasso

parfor i = 1:nt
    Gseg = GDat(tao(i):tao(i+1)-1,:);
    Kg(:,:,i) = BINCO_obs(Gseg,-5:0.1:5,100);
    Kg(:,:,i) = gaussIPF(Kg(:,:,i),(Gseg'*Gseg)/size(Gseg,1));
end