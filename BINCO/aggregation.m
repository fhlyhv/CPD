function [ntrue,copt]=aggregation(K_stb,nid)


%% predine
ntrue=0;
copt =0;
% p = size(K_stb,1);

%% compute v1 v2
Ku=triu(K_stb,1);
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
                if h2-h1>1e-5 && isnan(h1*h2)==0 %rule 3
                %% fit [v1,v2] to binomial distribution by minimizing KL divergence
                stbfit=stb(stb>=v1&stb<=v2);
                epdffit=epdf(stb>=v1&stb<=v2);
                if length(stbfit)>=4
                [a,b,pai] = bbinofit(stbfit,epdffit,nid);
                if b>1 %&& g>=1             
                    %% select threshold c for current lambda
                    c_array = flipud(stb(stb>=v2)); 
                    fdr_array = cumsum((1-pai)*bbinopdf(round(nid*c_array),nid,a,b))./cumsum(flipud(epdf(stb>=v2)));
%                     if any(fdr_array <= max(2/p/p-1,1e-2))
%                         idopt = find(fdr_array <= max(2/p/p-1,1e-2),1,'last');
%                         copt1 = c_array(idopt);
%                         copt2 = min(c_array(fdr_array==min(fdr_array)));
%                         copt = min(copt1,copt2);
%                     else
                        copt = min(c_array(fdr_array==min(fdr_array)));
%                     end
                    
                    %% compute number of true edges
                    nall=sum(Kc>=copt);
                    ntrue=nall*(1-min(fdr_array));
                end
                end
                end
            end
        end
    end
    end
end
end

    