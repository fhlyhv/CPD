1.To estimate the structure of the final graph, please call the function CGHVGM:

XDat=xlsread('nG(20)_2800');
K_est=CGHVGM(XDat);


Note that XDat should be a n x p dataset with n samples and p dimensions and K_est is the adjacent matrix corresponding to the resulting graphical model: K_est(i,j)=1 means there exists an edge between node i and node j, otherwise there is no edge between these two nodes.


2.If you wish to infer the number of hidden variables, please input

XDat=xlsread('nG(20)_2800');
[K_est,r]=CGHVGM(XDat,1);

The output r is the number of hidden variables. Besides, K_est here is the precision matrix corresponding to the graph. Compared with the true precision matrix (precision matrix20(169).xls), our result is a good estimation.  