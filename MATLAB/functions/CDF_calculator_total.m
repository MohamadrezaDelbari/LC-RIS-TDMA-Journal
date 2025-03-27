function [jj,ff,ff_dd] = CDF_calculator_total(f,Accuracy,MIN,MAX)
ff(size(f,2),1:Accuracy)=0;
for j=MIN:(MAX-MIN)/Accuracy:MAX
for k=1:size(f,2)  
for i=1:length(f)
    if  f(i,k)<j
        ff(k,round((j-MIN)*Accuracy/(MAX-MIN)))=ff(k,round((j-MIN)*Accuracy/(MAX-MIN)))+1;
    end
end
end
end
ff(1:size(f,2),:)=ff(1:size(f,2),:)./ff(1:size(f,2),length(ff));
jj=MIN+(MAX-MIN)/Accuracy:(MAX-MIN)/Accuracy:MAX;

for k=1:size(f,2) 
    ff_d(:,k) = diff(ff(k,:))./diff(jj);
end
 ff_dd=[zeros(1,k);ff_d]; 
end