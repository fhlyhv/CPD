function Kest=BINCO_obs(GDat,l,nuse)

% Bootstap Inference for Network COnstruction
% GDat is a nxp Gaussian distributed dataset, l is the initial region of
% log(lambda)
% Kest is the structure of precision matrix

% Bootstap Inference for Network COnstruction
% GDat is a nxp Gaussian distributed dataset, l is the initial region of
% log(lambda)
% Kest is the structure of precision matrix

p= size(GDat,2);
id= bootstrap_obs(GDat,nuse);
nid=size(id,2);
nl= length(l);
S=zeros(p,p,nid);
K_stb=zeros(p,p,nl);
ntrue_array=zeros(1,nl);
copt_array=zeros(1,nl);

%% compute stability for each edge for each lambda

for i=1:nid
    M=GDat(id(:,i),:);
    M(isnan(M))=0;
    S(:,:,i)=cov_normalize(cov(M,1));
end

clear GDat;


for i=1:nl
    K_stb(:,:,i)=stability(S,exp(l(i)),nid);
end



clear S;

K_stb(K_stb==0)=2;


%% find a range of lambda 

lb=1;
for i=1:nl
    Ku=triu(K_stb(:,:,i),1);
    Kc=Ku(:);
    Kc=Kc(Kc~=0);
    Kc(Kc==2)=0;
    [stb,epdf]=emppdf(Kc);
    ns = length(stb);
    if ns>3
        apdfint = mean(epdf(1:round(0.1*ns)));
        apdfend = mean(epdf(round(0.9*ns):ns));
        apdfmean = mean(epdf(round(0.4*ns):round(0.6*ns)));
        if apdfint-apdfmean>=2/length(Kc) && apdfend-apdfmean>=2/length(Kc);
        fspline=csaps(stb',epdf',0.999);
        spdf=fnval(fspline,0:1/nid:1);
        [~,id2]=min(spdf);
        v2=(id2-1)/nid;
        if v2<=0.8 && v2>min(stb) % rule 1
            [~,id1]=max(epdf(stb<=v2));
            v1=stb(id1(1));
            mu=(v1+v2)/2;
            h1=sum(epdf(stb>=v1&stb<=mu))/length(epdf(stb>=v1&stb<=mu));
            h2=sum(epdf(stb>=mu&stb<=v2))/length(epdf(stb>=mu&stb<=v2));
            if h1-h2>1e-5
                if epdf(id1(1))>epdf(abs(stb-v2-0.1/nid)==min(abs(stb-v2-0.1/nid))) && v1<v2 % rule 2
                    mu=(1+v2)/2;
                    h1=sum(epdf(stb>=v2&stb<=mu))/length(epdf(stb>=v2&stb<=mu));
                    h2=sum(epdf(stb>=mu))/length(epdf(stb>=mu));
                    if h2-h1>1e-5 && isnan(h1*h2)==0
                        lb=i;
                        break;
                    end
                end
            end
        end        
        end
    end
end

ub=nl;
for i=nl:-1:1
    Ku=triu(K_stb(:,:,i),1);
    Kc=Ku(:);
    Kc=Kc(Kc~=0);
    Kc(Kc==2)=0;
    [stb,epdf]=emppdf(Kc);
    ns = length(stb);
    if ns>3
        apdfint = mean(epdf(1:round(0.1*ns)));
        apdfend = mean(epdf(round(0.9*ns):ns));
        apdfmean = mean(epdf(round(0.4*ns):round(0.6*ns)));
        if apdfint-apdfmean>=2/length(Kc) && apdfend-apdfmean>=2/length(Kc);
        fspline=csaps(stb',epdf',0.999);
        spdf=fnval(fspline,0:1/nid:1);
        [~,id2]=min(spdf);
        v2=(id2-1)/nid;
        if v2<=0.8 && v2>min(stb) % rule 1
            [~,id1]=max(epdf(stb<=v2));
            v1=stb(id1(1));
            mu=(v1+v2)/2;
            h1=sum(epdf(stb>=v1&stb<=mu))/length(epdf(stb>=v1&stb<=mu));
            h2=sum(epdf(stb>=mu&stb<=v2))/length(epdf(stb>=mu&stb<=v2));
            if h1-h2>1e-5
                if epdf(id1(1))>epdf(abs(stb-v2-0.1/nid)==min(abs(stb-v2-0.1/nid))) && v1<v2 % rule 2
                    mu=(1+v2)/2;
                    h1=sum(epdf(stb>=v2&stb<=mu))/length(epdf(stb>=v2&stb<=mu));
                    h2=sum(epdf(stb>=mu))/length(epdf(stb>=mu));
                    if h2-h1>1e-5 && isnan(h1*h2)==0
                        ub=i;
                        break;
                    end
                end
            end
        end
        end
    end
end


%% compute maximum number of true edges for each lambda in the range

for i=lb:ub 
    [ntrue_array(i),copt_array(i)]=aggregation(K_stb(:,:,i),nid);
end

%% compute optiaml Kest    

[~,liopt]=max(ntrue_array);
copt=copt_array(liopt);
Kest=K_stb(:,:,liopt);
Kest(Kest==2|Kest<copt)=0;




