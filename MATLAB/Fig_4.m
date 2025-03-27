%%----------------------------------------------------%%
%%----- Mohamadreza Delbari
%%----- Please cite: Fast Reconfiguration of LC-RISs in TDMA Protocol
%%----- Please cite: DOI: 10.1109/ICCWorkshops59551.2024.10615422
%%----------------------------------------------------%%
clear all
clc
close all
t=0:160;
f=t*0;
for t=10:60
    f(t+1)=360*(1-exp(-(t-10)/9));
end
for t=60:70
    f(t+1)=360;
end
for t=70:160
    f(t+1)=360*(exp(-(t-70)/29));
end
t=0:160;
plot(t,f)
AA=[t',f'];
