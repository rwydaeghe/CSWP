fas1=zeros(n_of_samples,1);
maxf=1/dt;
df=maxf/n_of_samples;
for loper=1:n_of_samples
   fas1(loper,1)=df*(loper-1);
end
k=2*pi*fas1/c
f=zeros(n_of_samples,1);
for i=1:n_of_samples
    f(i)=analytic(1,1,k(i),10,2,0);
end
f
analytic(1,1,k(8192),10,2,0)
    

