function [SNR,tt,tsnr,LCT] = SNR_vs_time(Param,omega,xmax,xmin)

p_mu=Param.p_mu;
y_step=Param.y_step;
K=length(omega(1,:));
LCT(1:K)=10^10;
threshold(1:K)=Param.THR;
tsnr=0;
N=length(omega);

tau_neg = 29; % bottleneck
tau_pos = 9; % positive
omega_1=omega(:,K);
omega_2=omega(:,1);
omega_temp=zeros(N,1);
II=20;
JJ=K;
SNR((JJ-1)*II+II)=0;
tt((JJ-1)*II+II)=0;
t(II)=0;
for jj=1:JJ
    for ii=1:II
        t(ii)=(ii-1)*3;
        parfor i=1:N
            if omega_1(i)<omega_2(i)
                omega_temp(i)=omega_1(i)+(xmax-omega_1(i))*(1-exp(-t(ii)/(tau_pos)));
                if omega_temp(i)>omega_2(i)
                    omega_temp(i)=omega_2(i);
                end
            else
                omega_temp(i)=omega_1(i)+(xmin-omega_1(i))*(1-exp(-t(ii)/(tau_neg)));
                if omega_temp(i)<omega_2(i)
                    omega_temp(i)=omega_2(i);
                end
            end
        end
        Gamma=diag(exp(1i*omega_temp));
        %H_eff = H_d+H_r*Gamma*H_i;
        for i=1:length(Param)
            Param(i).y_start=p_mu(jj,2)-5*y_step;
            Param(i).y_end=p_mu(jj,2)+5*y_step;
        end
        SNR((jj-1)*II+ii) = min(SNR_calculation(Param,Gamma,0));
        if ii>1 && ii<length(SNR)-1
            if (SNR((jj-1)*II+ii-1)<threshold(jj)) && (SNR((jj-1)*II+ii)>threshold(jj))
                %LCT(jj)=(II-1)/(II/2)*(jj-1)+t(ii);
                LCT(jj)=t(ii-1)+(10-SNR((jj-1)*II+ii-1))/(SNR((jj-1)*II+ii)-SNR((jj-1)*II+ii-1))*3;
            end
        end

        tt((jj-1)*II+ii)=(jj-1)*(II-1)*3+t(ii);

    end
    if jj<K
        omega_1=omega(:,jj);
        omega_2=omega(:,jj+1);
    end
    if (LCT(jj)>10^9) && (SNR((jj-1)*II+ii)>threshold(jj))
        LCT(jj)=0;
    end

    %[omega_1, omega_2]=deal(omega_temp,mod(jj,2)*W1+(1-mod(jj,2))*W2);
end



end