
tao = unique([1,tao,size(GDat,1)]);
nt = length(tao)-1;
Kh = zeros(p,p,nt); % hidden variable
nh = zeros(1,nt);
parfor i = 1:nt
    Gseg = GDat(tao(i):tao(i+1)-1,:);
    [Kh(:,:,i),nh(i)] = GHVGM(Gseg,1);
end
