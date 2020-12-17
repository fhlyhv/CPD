function [a,b,pai] = bbinofit(stbfit, epdffit, nuse)

stbfitn = round(stbfit*nuse);

opts = optimset('Display','off','TolX',1e-6,'TolFun',1e-6);
p = lsqnonlin(@KL_divergence,[0.5;2;0.1],[0;1;0],[1;1e3;1],opts);

a = p(1);
b = p(2);
pai = p(3);
% pai = 1-2/sum(bbinopdf([0:nuse-1,1:nuse]',nuse,p(1),p(2)));

% --------------------------------------------------------------------
function kldiv = KL_divergence(p)

kldiv = epdffit-(1-p(3))*bbinopdf(stbfitn,nuse,p(1),p(2));




end

end