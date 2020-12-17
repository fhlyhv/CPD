function Id = bootstrap_obs(XDat,N)




% n=size(XDat,1);
% Id=zeros(round(n/2),N);
% 
% for i=1:N 
%     I=floor(rand(n,1)*n)+1;
%     I(I==n+1)=n;
%     Id(:,i)=I;
% end
% 
% % for i = 1:2:N-1
% %     I = randperm(n);
% %     Id(:,i) = I(1:round(n/2));
% %     Id(:,i+1) = I(round(n/2)+1:end);
% % end


n=size(XDat,1);
hn = floor(n/2);
Id=zeros(hn,N);

for i=1:2:N-1
   I = randperm(n);
   Id(:,i) = I(1:hn);
   Id(:,i+1) = I(hn+1:2*hn);
%     I=floor(rand(n,1)*n)+1;
%     I(I==n+1)=n;
%     Id(:,i)=I;
end

