%% generate aritficial data
p = 25;
T = 1900; %2400;
ns = 4;


cptrue = [400+randperm(100,1),900+randperm(100,1),1400+randperm(100,1)]; % true change point position
%[200+randperm(300,1),700+randperm(300,1),1200+randperm(300,1),1700+randperm(300,1)];
fprintf('The true change points are: [ ');
fprintf('%g ', cptrue);
fprintf(']\n');

cp2 = [1,cptrue;cptrue-1,T];
GDat = zeros(T,p);
Ktrue = zeros(p,p,ns);
for j = 1:ns
    [GDat(cp2(1,j):cp2(2,j),:),Ktrue(:,:,j)] = ArtiDatGen(p,cp2(2,j)-cp2(1,j)+1,0.05); %0.05+0.15*rand(1)
end
XDat = G2NG(GDat);

%% Detect change point

GDat = norminv(empcdf_con(XDat)); % convert to Gaussian data
cpca = APELT(GDat,100,200); % estimated change points
fprintf('The estimated change points are: [ ');
fprintf('%g ', cpca);
fprintf(']\n');

%% learn graphical model structure
addpath(genpath('.\BINCO'));
tau = unique([1,cpca,T]);
nt = length(tau)-1;
KBINCO = zeros(p,p,nt);
for j = 1:nt
    Gseg = GDat(tau(j):tau(j+1)-1,:);
    KBINCO(:,:,j) = BINCO_obs(Gseg,-5:0.1:4,100);
end

%% check graph estimation accuracy
prcs = (sum(Ktrue(:)~=0 & KBINCO(:)~=0)-ns*p)/(sum(KBINCO(:)~=0)-ns*p);    
rc = (sum(Ktrue(:)~=0 & KBINCO(:)~=0)-ns*p)/(sum(Ktrue(:)~=0)-ns*p);
f1s = 2*prcs*rc/(prcs+rc);

fprintf('precision = %d, recall = %d, f1-score = %d\n', prcs, rc, f1s); 