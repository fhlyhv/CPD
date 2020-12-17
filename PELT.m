function [tao,Fn] = PELT(GDat,ndefault,b)

% PELT: pruned exact linear time method
% GDat - Gaussian distributed data
% ndefault - minimum number of time points in each segment
% b - penalty parameter
% Yu Hang, NTU, Jun, 2013



% GDat = GDat(1:n,:);
n = size(GDat,1);
cp = zeros(1,n);

% F0 = -b;
F = zeros(1,n);
R = 0:n-1;
Rflag = [true,false(1,ndefault-1),true(1,n-ndefault)];

for i = ndefault:n
    Rtemp = R(Rflag(1:i-ndefault+1));  %(1:nR)
    nR = length(Rtemp);
    Fr = zeros(1,nR);
    for j = 1:nR
        if Rtemp(j) == 0
            Fr(j) = Contra(GDat(1:i,:),n); %F0+Contra(GDat(1:i,:),n)+b;
        else
            Fr(j) = F(Rtemp(j))+Contra(GDat(Rtemp(j)+1:i,:),n)+b;
        end
    end
    [F(i),id] = min(Fr);  %(1:nR)
    cp(i) = Rtemp(id);
    %if i <n
    %   Rflag(i+1)=1;
    Rflag(Rtemp(Fr-b>=F(i))+1) = false;  %(1:nR)
    %end
end

tao = [];

while cp(i)~=0
    tao = cat(2,cp(i),tao);
    i = cp(i);
end

% tao = tao(tao>0);

Fn = F(n);
if isempty(tao)
    tao = 0;
    Fn = -Inf;
end


