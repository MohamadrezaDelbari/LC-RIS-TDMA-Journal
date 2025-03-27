%%----------------------------------------------------%%
%%----- Mohamadreza Delbari
%%----- Please cite: Fast Reconfiguration of LC-RISs in TDMA Protocol
%%----- Please cite: DOI: 10.1109/ICCWorkshops59551.2024.10615422
%%----------------------------------------------------%%

function [gamma,deltagamma,avgt] = algorithm(Param,lambda,gamma,alpha,Imax,SNR_thr,tmax,xmax,xmin)
p_mu=Param.p_mu;
y_step=Param.y_step;
K=length(gamma(1,:));
lambda(1:K)=lambda;
delta_tmax=tmax/Imax;
for ii=1:Imax
    for k=1:K
        if k==1
            omega(:,k) = fminfinder(gamma(:,K),gamma(:,k+1),  lambda(k), gamma(:,k),tmax,xmax,xmin);
            for i=1:length(Param)
                Param(i).y_start=p_mu(k,2)-5*y_step;
                Param(i).y_end=p_mu(k,2)+5*y_step;
            end
            SNR=SNR_calculation(Param,diag(exp(1i*omega(:,k))),0);
            if min(SNR)<SNR_thr
                lambda(k)=lambda(k)/alpha;
            else
                lambda(k)=lambda(k)*alpha;
                gamma(:,k)=omega(:,k);
            end
        elseif k==K
          omega(:,k) = fminfinder(gamma(:,K-1),gamma(:,1),  lambda(k), gamma(:,k),tmax,xmax,xmin);
            for i=1:length(Param)
                Param(i).y_start=p_mu(k,2)-5*y_step;
                Param(i).y_end=p_mu(k,2)+5*y_step;
            end
            SNR=SNR_calculation(Param,diag(exp(1i*omega(:,k))),0);
            if min(SNR)<SNR_thr
                lambda(k)=lambda(k)/alpha;
            else
                lambda(k)=lambda(k)*alpha;
                gamma(:,k)=omega(:,k);
            end
        else
            omega(:,k) = fminfinder(gamma(:,k-1), gamma(:,k+1) , lambda(k), gamma(:,k),tmax,xmax,xmin);
            for i=1:length(Param)
                Param(i).y_start=p_mu(k,2)-5*y_step;
                Param(i).y_end=p_mu(k,2)+5*y_step;
            end
            SNR=SNR_calculation(Param,diag(exp(1i*omega(:,k))),0);
            if min(SNR)<SNR_thr
                lambda(k)=lambda(k)/alpha;
            else
                lambda(k)=lambda(k)*alpha;
                gamma(:,k)=omega(:,k);
            end
        end
    end

    % convergence rate
    deltagamma(ii,1)=(norm(gamma(:,1)-gamma(:,2),2)^2+norm(gamma(:,1)-gamma(:,3),2)^2+norm(gamma(:,2)-gamma(:,3),2)^2)/3/length(gamma);
    maxgamma(ii,1)=max([max(abs(gamma(:,1)-gamma(:,2)));max(abs(gamma(:,1)-gamma(:,3)));max(abs(gamma(:,2)-gamma(:,3)))]);
    for k=1:K
        if k==1
            t(ii,k)=max(9*(sign(gamma(:,k)-gamma(:,K))+1)/2.*log((xmax-gamma(:,K))./(xmax-gamma(:,k)))-29*(sign(gamma(:,k)-gamma(:,K))-1)/2.*log((xmin-gamma(:,K))./(xmin-gamma(:,k))));
        else
            t(ii,k)=max(9*(sign(gamma(:,k)-gamma(:,k-1))+1)/2.*log((xmax-gamma(:,k-1))./(xmax-gamma(:,k)))-29*(sign(gamma(:,k)-gamma(:,k-1))-1)/2.*log((xmin-gamma(:,k-1))./(xmin-gamma(:,k))));
        end
    end
    tmax=tmax-delta_tmax;
end
 avgt=sum(t')'/K;
end
