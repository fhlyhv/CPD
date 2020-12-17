%% generate aritficial data
p = 25;
T = 1900; %2400;
ns = 4;

cptrue = [400+randperm(100,1),900+randperm(100,1),1400+randperm(100,1)];
%[200+randperm(300,1),700+randperm(300,1),1200+randperm(300,1),1700+randperm(300,1)];
fprintf('The true change points are: [ ');
fprintf('%g ', cptrue);
fprintf(']\n');

cp2 = [1,cptrue;cptrue-1,T];
GDat = zeros(T,p);
Ktrue = zeros(p,p,ns);
nh = zeros(ns,1);
for j = 1:ns
    [~,K] = ArtiDatGen(p,1,0.05); %0.05+0.15*rand(1)
    Ktrue(:,:,j) = K;
    nh(j) = randperm(2,1);
    for k = 0:nh(j)-1
        Koh = zeros(p+k,1);
        nho = ceil((p+k)*(0.8+0.2*rand(1)));
        Koh(randperm(p+k,nho)) = sign(rand(nho,1)-0.5).*(0.5*rand(nho,1));
        Khh = abs(sum(Koh))+4;
        K = [K,Koh;Koh',Khh];
    end
    S = inv(K);
    S = S(1:p,1:p);
    GDat(cp2(1,j):cp2(2,j),:) = mvnrnd(zeros(p,1),S,cp2(2,j)-cp2(1,j)+1);
end
XDat = G2NG(GDat);


%% Detect change point

GDat = norminv(empcdf_con(XDat)); % convert to Gaussian data
cpca = APELT(GDat,100,200); % estimated change points
fprintf('The estimated change points are: [ ');
fprintf('%g ', cpca);
fprintf(']\n');

%% learn graphical model structure
addpath(genpath('.\CGHVGM with ss'));
tau = unique([1,cpca,T]);
nt = length(tau)-1;
Kss = zeros(p,p,nt); % conditional precision matrices
nhss = zeros(1,nt); % no. of hidden variables

for j = 1:nt
    Gseg = GDat(tau(j):tau(j+1)-1,:);
    [Kss(:,:,j),nhss(j)] = GHVGM(Gseg,1);
end

%% check graph estimation accuracy
prcs = (sum(Ktrue(:)~=0 & Kss(:)~=0)-ns*p)/(sum(Kss(:)~=0)-ns*p);    
rc = (sum(Ktrue(:)~=0 & Kss(:)~=0)-ns*p)/(sum(Ktrue(:)~=0)-ns*p);
f1s = 2*prcs*rc/(prcs+rc);

fprintf('precision = %d, recall = %d, f1-score = %d\n', prcs, rc, f1s); 