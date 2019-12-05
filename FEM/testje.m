max=30;
data=zeros(1,max);
for N=3:15
    N
    %figure('Name',sprintf('%d', N))
    tic
    ComputeEigsTE(N);
    data(N)=toc;
end
%hold on
%scatter(1:max-1,diff(sqrt(data)))