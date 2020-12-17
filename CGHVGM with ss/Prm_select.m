function [X1,X2]=Prm_select(Sigma,K_est,lambda)


%%*************************************************************************
%% Solve: 
%% min { <X1,Sigma> - logdetX1 + beta*<I,X2> + M<K_est,Xs>
%%
%% X1 + X2 - Xs = 0. 
%% X1 positive definite, X2 positive semidefinite, Xp,Xm >= 0.
%%
%%*************************************************************************

%    HOME = '.\LogdetPPA-0'; 
%    addpath(strcat(HOME,'/solver/'))
%    addpath(strcat(HOME,'/solver/mexfun'))
%    addpath(strcat(HOME,'/util/'))
%    ttime  = clock;
%% set up SDP data in SDPT3 format
%%
      n=size(Sigma,1);
      invD = speye(n,n); 
      n2 = n*(n+1)/2; 
      Identity = speye(n2); 
      b = zeros(n2,1);   
      C{1} = Sigma; 
      blk{1,1} = 's'; blk{1,2} = n;
      
      At{1} = Identity; 
      %% 
      blk{2,1} = 's'; blk{2,2} = n;
      At{2,1}  = Identity; 
      C{2,1}   = lambda*speye(n,n); 
      %%
      blk{3,1} = 'l'; blk{3,2} = 2*n2;  
      At{3,1} = [-Identity,Identity]';
      [row, col] = find(triu(K_est==1)); 
      idx = row + col.*(col-1)/2;
      ee  = 1e4*ones(n2,1);
      ee(idx) = zeros(length(idx),1); 
      C{3,1} = [ee;ee];
%       fprintf('\n Set up data time = %3.2f',etime(clock,ttime)); 
      runPPA = 1; 
      if (runPPA)
         OPTIONS.smoothing  = 1;
         OPTIONS.scale_data = 0; %% or 2;
         OPTIONS.plotyes    = 0; 
         OPTIONS.tol        = 1e-6;
         mu = [1; 0; 0];
         [obj,X,y,Z,runhist] = logdetPPA(blk,At,C,b,mu,OPTIONS);
         %obj = sum(sum(Sigma.*X{1}))-sum(log(eig(X{1})))+rho*sum(sum(abs(X{1}))); 
         X1 = invD*X{1}*invD; X1 = 0.5*(X1+X1');   
         X2 = invD*X{2}*invD; X2 = 0.5*(X2+X2');   
      end
%%*************************************************************************
