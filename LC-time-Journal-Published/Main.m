%%----------------------------------------------------%%
%%----- Mohamadreza Delbari
%%----- Please cite: Fast Reconfiguration of LC-RISs in TDMA Protocol
%%----- Please cite: DOI: 10.1109/ICCWorkshops59551.2024.10615422
%%----------------------------------------------------%%

% Delete all the path cashes
restoredefaultpath

%% Parameters
close all
clear variables
clc
addpath('./functions')
load('Data/workspace.mat')
%% channel
ParamC = struct('p_bs',p_bs,'p_irs',p_irs,'p_mu',p_mu,...
            'factor_Hd',factor_Hd,'factor_Hi',factor_Hi,'factor_Hr',factor_Hr,...
            'eta',eta,'d0',d0,'beta',beta,'h_blk',h_blk,...
            'lambda',lambda,'N_bs_x',N_bs_x,'N_bs_z',N_bs_z,'N_mu',N_mu,...
            'N_y',N_y,'N_z',N_z,...
            'd_bs_x',d_bs_x,'d_bs_z',d_bs_z,'d_irs_y',d_irs_y,'d_irs_z',d_irs_z,...
            'Krice_dB',Krice_dB,'Vscatter',Vscatter,'f',f,'S',S);

        [H_d,H_i,H_r,Param_output] = func_channel(ParamC);


%% Initial phase shifter design according to arxiv.org/pdf/2401.08237
for kk=1:K
    if kk==1
        pp_mu{kk}=[p_mu(kk,:)-[0 1 0];p_mu(kk,:)+[0 1 0]];
    elseif kk==2
        pp_mu{kk}=[p_mu(kk,:)-[0 1 0];p_mu(kk,:)+[0 1 0]];
    elseif kk==3
        pp_mu{kk}=[p_mu(kk,:)-[0 1 0];p_mu(kk,:)+[0 1 0]];
    end
    pris=Param_output.pp_ris_g;
    [WW{kk}] = NF_RIS_design(f,N_y,N_z,pp_mu{kk},p_bs,pris);
end
%save("W.mat","WW")
%WW=load("W.mat");
%WW=WW.WW;

reverse=1; % Impact of ordering in user serving
if reverse==1
    p_mu_temp=p_mu(2,:);
    p_mu(2,:)=p_mu(3,:);
    p_mu(3,:)=p_mu_temp;
    WW_temp=WW{2};
    WW{2}=WW{3};
    WW{3}=WW_temp;
end

%% Test receive SNR
Param=struct('p_bs',p_bs,'p_irs',p_irs,'p_mu',p_mu,...
            'factor_Hd',factor_Hd,'factor_Hi',factor_Hi,'factor_Hr',factor_Hr,...
            'eta',eta,'d0',d0,'beta',beta,'h_blk',h_blk,...
            'lambda',lambda,'N_bs_x',N_bs_x,'N_bs_z',N_bs_z,'N_mu',N_mu,...
            'N_y',N_y,'N_z',N_z,...
            'd_bs_x',d_bs_x,'d_bs_z',d_bs_z,'d_irs_y',d_irs_y,'d_irs_z',d_irs_z,...
            'Krice_dB',Krice_dB,'Vscatter',Vscatter,'f',f,'S',S,'x_start',p_mu(1,1),'x_end',p_mu(1,1),'y_start',-10,'y_end',10,...
            'x_step',0.2,'y_step',0.2,'K',K,'Power',Power,'THR',10);

for kk=1:K
SNR(:,:,kk)=SNR_calculation(Param,WW{kk},0);
end

%% Initial Phase shift

for kk=1:K
    omega_init(:,kk)=diag(WW{kk});
end
gamma=angle(omega_init)+pi;

%% New design
xmin=0;
xmax=4*pi;
Imax=35;
lambda=[1,1,1];
tmax=100;
THR=10.3;
[gamma_new,deltagamma,avgt]=algorithm(Param,lambda,gamma,0.95,Imax,THR,tmax,xmax,xmin);
AAA_temp=[(1:length(deltagamma))',deltagamma];

for kk=1:K
   WW{kk}=diag(exp(1i*gamma_new(:,kk)));
   SNR(:,:,kk)=SNR_calculation(Param,WW{kk},0);
end

%% Generating Fig. 6
    i=1:Imax;
    figure
    plot(i,avgt)
    AAA_6=[i',avgt];
%% Generating Fig. 7
for kk=1:K
    if kk==1
        gamma_delta(:,kk)=gamma(:,kk)-gamma(:,K);
        gamma_new_delta(:,kk)=gamma_new(:,kk)-gamma_new(:,K); 
    else
        gamma_delta(:,kk)=gamma(:,kk)-gamma(:,kk-1);
        gamma_new_delta(:,kk)=gamma_new(:,kk)-gamma_new(:,kk-1); 
    end
end
[gamma_index,gamma_cdf,gamma_pdf]=CDF_calculator_total(gamma_delta,50,-2*pi,2*pi);
[gamma_new_index,gamma_new_cdf,gamma_new_pdf]=CDF_calculator_total(gamma_new_delta,50,-2*pi,2*pi);

figure
plot(gamma_index,gamma_pdf,gamma_new_index,gamma_new_pdf)
AAA_7=[gamma_index',gamma_pdf,gamma_new_pdf];
AAA_7=[[-2*pi zeros(1,6)]',AAA_7']';

%% Generating Fig. 8
Param(1).THR=10;

[SNRR,t,tsnr,LCT] = SNR_vs_time(Param,gamma,xmax,xmin);
[SNRR_new,t,tsnr,LCT_new] = SNR_vs_time(Param,gamma_new,xmax,xmin);
figure
plot(t,SNRR,t,SNRR_new)
THR=ones(1,length(t))*10;
AAA_8=[t',SNRR',SNRR_new',THR'];

%% Generating Fig .9
LCT_mean=mean(LCT);
LCT_new_mean=mean(LCT_new);

Ts_e=(0:0.1:3);
Ts=10.^Ts_e;
data_rate=max((Ts-LCT_mean)./Ts*log2(1+10),0);
data_rate_new=max((Ts-LCT_new_mean)./Ts*log2(1+10),0);
figure
semilogx(Ts,data_rate,Ts,data_rate_new)
AAA_9=[Ts' data_rate' data_rate_new'];