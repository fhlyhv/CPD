function [XDat,Ktrue] = ArtiDatGen(P,N,r)

% N -- sample size
% P -- dimension
% r -- no. of non-zero elements in the Cholesky decomposition of the
% precision matrix



[rid,cid] = find(triu(ones(P),1));
Pe = length(rid);

Pt = ceil(r*Pe); 
ide = randperm(Pe,Pt).';

Ku = sparse(rid(ide),cid(ide),sign(rand(Pt,1)-0.5).*(0.5+0.5*rand(Pt,1)),P,P);
Ku = spdiags(1+0.5*rand(P,1),0,Ku);
Su = speye(P)/Ku;
Sd = spdiags(sqrt(sum(Su.^2,2)),0,P,P);
Ku = Ku*Sd;
Ktrue = Ku.'*Ku;
XDat = randn(N,P)/Ku.';