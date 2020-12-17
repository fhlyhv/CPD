function [K_est,r] = GHVGM(GDat,flag)

% XDat should be a n*p datset, with n samples and p dimensions.

%% Gaussian copula
% GDat=copula(XDat);
[n,p] = size(GDat);
S=GDat'*GDat/n; %cov(GDat,1);
EV = p*(p-1)/2*0.1;

%% structure learning
[l_r,g_r,K_array]=rangsel(GDat,0.01,0.01,0.3,0.1,0.1,3);
K_est=stasel_rnran(GDat,100,0.2,0.2,l_r(1),0.01,l_r(2),g_r(1),0.1,g_r(2),EV);

%% parameter learning 
if nargin == 2 && flag ==1
    lambda = lambda_sel(K_est,K_array,0.01,0.01,0.3,0.1,0.1,3);
    [X_1,X_2]=Prm_select(S,K_est,lambda);
    K=X_1+X_2;
    K(abs(K)<1e-3)=0;
    K_est=K;
    r=rank(X_2);
end

