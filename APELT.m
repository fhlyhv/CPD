function [tao,taoall] = APELT(GDat,ndefault,nb)

% Adaptive Pruned Exact Linear-Time Method
% tao -- change point position
% taoall -- change point position for all the selected beta
% GDat -- Gaussian distributed data
% ndefault -- default minimum number of samples in each time segment
% nb -- number of testing beta
% Yu Hang, Kyoto University, Nov, 2012

%% predine
[n,p] = size(GDat);
taoall = cell(1,nb);
F = zeros(1,nb);
ntaoall = zeros(1,nb);
tb = 0.5*p*(p+1)*log(n)/n/4; %(1+(1:nb)*0.02)*p^2*log(n)/n;  %0.01
tao = [];

while 1 
    [taoall{1},F(1)] = PELT(GDat,ndefault,tb);
    ntaoall(1) = length(taoall{1});
    if ntaoall(1)>10
        tb = tb*1.5;
    elseif ntaoall <5
        tb = tb*0.8;
    else
        tb = tb+(1:nb-1)*0.1*p*(p+1)*log(n)/n/4;
        break;
    end
end

%% identify change points for different values of penalty parameter
for i=2:nb
    [taoall{i},F(i)] = PELT(GDat,ndefault,tb(i));
    ntaoall(i) = length(taoall{i});
    tbl = tabulate(ntaoall);
    id = find(tbl(:,1)==ntaoall(i)+1);
    if ~isempty(id)
        tmp = tbl(id,2);
    else
        tmp = 0;
    end
    if tbl(tbl(:,1)==ntaoall(i),2)-tmp>=2 && ntaoall(i)<=10 
        tao = taoall{i};
        break;
    end
end

% %% find the best beta
% nt = max(ntaoall):-1:min(ntaoall);
% nnt = length(nt);
% lt = zeros(1,nnt);
% for i = 1:nnt
%     lt(i) = sum(ntaoall == nt(i));
% end
% dlt = lt(3:end)-lt(2:end-1);
% id = find(dlt>2,1,'first');
% 
% nc = nt(id+2);
% idt = find(ntaoall == nc);
% tao = taoall{idt(1)};


