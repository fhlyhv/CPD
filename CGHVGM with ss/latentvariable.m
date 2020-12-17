function [X1,X2]=latentvariable(data,lambda,gamma,d)


%%*************************************************************************
%% Solve: 
%% min { <X1,Sigma> - logdetX1 + beta*<I,X2> + rho*<E,Xp> + rho*<E,Xm> 
%%
%% X1 + X2 - Xp + Xm = 0. 
%% X1 positive definite, X2 positive semidefinite, Xp,Xm >= 0.
%%
%%*************************************************************************

%    HOME = '.\LogdetPPA-0'; 
%    addpath(strcat(HOME,'/solver/'))
%    addpath(strcat(HOME,'/solver/mexfun'))
%    addpath(strcat(HOME,'/util/'))
%    ttime  = clock;
%%
%% Data setting up
%%
   rho = lambda*gamma;
   beta = lambda;
   Gene = data; 

   Tmp = isnan(Gene);
   Gene1 = Gene;
   for i = 1:size(Gene,1)
      for j = 1:size(Gene,2)
         if (Tmp(i,j)==1); Gene1(i,j)=0; end;
      end
   end
%%
   dim = [d,500,1000,1500]; 
   for kk = [1:1]; 
      n = dim(kk); 
     % vv = var(Gene1'); 
     %[vvsort,idx] = sort(vv); 
     % vvsort = vvsort(end:-1:1); idx = idx(end:-1:1); 
      %% select n variables with highest variances
     % index = idx(1:n); 
     % Sigma = cov(Gene1(index,:)');                                     
      Sigma=cov(Gene1',1);
      Sigma = 0.5 * (Sigma + Sigma');
%%
      Sigma_org = Sigma; invD = speye(n,n); 
      %%dd = max(1,sqrt(diag(Sigma))); invD = diag(1./dd); 
      %%Sigma = invD*Sigma*invD; 
      %%Sigma = 0.5*(Sigma+Sigma');
%%
%% set up SDP data in SDPT3 format
%%
      n2 = n*(n+1)/2; 
      b = zeros(n2,1);   
      C{1} = Sigma; 
      blk{1,1} = 's'; blk{1,2} = n;
      [Iall,Jall] = find(triu(ones(n)));
      tmp = [Iall,Jall]; 
      m2 = size(tmp,1); 
      Icomp = tmp(:,1); Jcomp = tmp(:,2);
      Itmp  = Icomp + Jcomp.*(Jcomp-1)/2; 
      Atmp  = spconvert([Itmp,[1:m2]',ones(m2,1); n2,m2,0]); 
      At{1} = Atmp; 
      %% 
      blk{2,1} = 's'; blk{2,2} = n;
      At{2,1}  = Atmp; 
      C{2,1}   = beta*speye(n,n); 
      %%
      blk{3,1} = 'l'; blk{3,2} = 2*n2;  
      Identity = speye(n2); 
      At{3,1} = [-Identity,Identity]';
      idx = find(Icomp == Jcomp); 
      ee  = sqrt(2)*ones(m2,1); 
      if ~isempty(idx); ee(idx) = ones(length(idx),1); end
      C{3,1} = rho*[ee; ee];
%       fprintf('\n Set up data time = %3.2f',etime(clock,ttime)); 
      runPPA = 1; 
      if (runPPA)
         OPTIONS.smoothing  = 1;
         OPTIONS.scale_data = 0; %% or 2;
         OPTIONS.plotyes    = 0; 
         OPTIONS.tol        = 1e-6;
         mu = [1; 0; 0];
         [obj,X,y,Z,runhist] = logdetPPA(blk,At,C,b,mu,OPTIONS);
         obj = sum(sum(Sigma.*X{1}))-sum(log(eig(X{1})))+rho*sum(sum(abs(X{1}))); 
         X1 = invD*X{1}*invD; X1 = 0.5*(X1+X1');   
         X2 = invD*X{2}*invD; X2 = 0.5*(X2+X2');   
      end
   end
%%*************************************************************************
