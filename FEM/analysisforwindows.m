data=zeros(7,2)
for i=1:30
    i
    if i==1|i==2
        data(i,:)=[0,0]
    else
        data(i,:)=ComputeEigsTE(i)
    end
end
plot([3:30],data(3:end,:))

%hiervoor moet je output=ComputeTE doen en uncommenten daarin