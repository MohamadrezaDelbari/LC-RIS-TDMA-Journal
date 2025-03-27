%% Edited by Mohamadreza
function [H_d,H_i,H_r,Param_output] = func_channel(Param)

p_bs = Param.p_bs; % BS's positions
p_ris = Param.p_irs; % RIS's positions
p_mu = Param.p_mu; % MU's positions

N_bs_x = Param.N_bs_x; % number of BS's antenna along x axis
N_bs_z = Param.N_bs_z; % number of BS's antenna along z axis
N_bs = N_bs_x*N_bs_z;

d_bs_x = Param.d_bs_x; % BS element spacing along x axis
d_bs_z = Param.d_bs_z; % BS element spacing along z axis

N_mu = Param.N_mu; % number of MU's antenna

N_y = Param.N_y; % number of IRS tiles along y axis
N_z = Param.N_z; % number of IRS tiles along z axis
N = N_y*N_z;

d_ris_y = Param.d_irs_y; % RIS element spacing along y axis
d_ris_z = Param.d_irs_z; % RIS element spacing along z axis

Krice_dB = Param.Krice_dB; % K factor 

factor_Hd = Param.factor_Hd; % normalization factor
factor_Hi = Param.factor_Hi; % normalization factor
factor_Hr = Param.factor_Hr; % normalization factor

eta = Param.eta; % pathloss exponent
d0 = Param.d0; % reference distance
beta = Param.beta; % reference loss
h_blk = Param.h_blk; % blockage gain for dircet link

lambda = Param.lambda; % wavelength
f=Param.f; % frequency

[K,~] = size(p_mu); % number of users

eta_d = eta(1); % pathloss exponent for BS-MU link
eta_i = eta(2); % pathloss exponent for BS-IRS link
eta_r = eta(3); % pathloss exponent for IRS-MU link

kappa = 2*pi/lambda; % wavenumber

Krice = db2pow(Krice_dB); % Rician factor
Krice_d = Krice(1); % Rician factor for BS-MU link
Krice_i = Krice(2); % Rician factorfor BS-IRS  link
Krice_r = Krice(3); % Rician factor for IRS-MU link

% non-LOS scatterers
Vscatter_d = Param(1).Vscatter; % scattering volume for paths in BS-MU link
Vscatter_i = Param(2).Vscatter; % scattering volume for paths in BS-MU link
Vscatter_r = Param(3).Vscatter; % scattering volume for paths in BS-MU link
S=Param.S;

%% large scale pathloss and fading:
for kk=1:K
    d = norm(p_bs-p_mu(kk,:));
    hbar_pl_d(kk) = beta*(d0/d)^eta_d; % pathloss

    d = norm(p_ris-p_mu(kk,:));
    hbar_pl_r(kk) = beta*(d0/d)^eta_r; % pathloss
end

d = norm(p_bs-p_ris);
hbar_pl_i = beta*(d0/d)^eta_i; % pathloss

%% Positions
% antenna and element positions:
ParamP = struct('lambda',lambda,'N_bs_x',N_bs_x,'N_bs_z',N_bs_z,'p_bs',p_bs,'p_ris',p_ris,...
        'N_y',N_y,'N_z',N_z,'d_bs_x',d_bs_x,'d_bs_z',d_bs_z,'d_ris_y',d_ris_y,'d_ris_z',d_ris_z);
[pp_bs,pp_ris]  = func_bs_irs_ant_p(ParamP);
pp_bs_g=p_bs+pp_bs;
pp_ris_g=p_ris+pp_ris;

%% Small scale fading

%---------------------------BS_RIS---------------------------
% scattering positions:
for s=1:S
    pp_scr_i(s,:) = (Vscatter_i(2,:)-Vscatter_i(1,:)).*rand(1,3)+Vscatter_i(1,:); % scattering cluster position
    for kk=1:K
        % scattering positions:
        pp_scr_d{kk}(s,:) = (Vscatter_d(2,:)-Vscatter_d(1,:)).*rand(1,3)+Vscatter_d(1,:); % scattering cluster position
        pp_scr_r{kk}(s,:) = (Vscatter_r(2,:)-Vscatter_r(1,:)).*rand(1,3)+Vscatter_r(1,:); % scattering cluster position
    end
end

%---------------------------RIS-Users---------------------------


%% Steering vectors:

% BS-IRS steering vector:
a_bs_ris = zeros(N_bs,1+S);
d_BStoRIS = (p_ris-p_bs)/norm(p_ris-p_bs);
a_bs_ris(:,1) = exp(1j*kappa*pp_bs*d_BStoRIS');
for s=1:S
    d_BStoScr = (pp_scr_i(s,:)-p_bs)/norm(pp_scr_i(s,:)-p_bs);
    a_bs_ris(:,s+1) = beta*exp(1j*kappa*pp_bs*d_BStoScr')/(norm(pp_scr_i(s,:)-p_bs)*norm(pp_scr_i(s,:)-p_ris))^eta_i;
end


% IRS-BS steering vectors:
a_ris_bs_vec = [];
    a_ris_bs(:,1) = exp(1j*mod(kappa*vecnorm(p_bs-pp_ris_g,2,2),2*pi));
for s=1:S
    a_ris_bs(:,1+s) = beta*exp(1j*mod(kappa*vecnorm(pp_scr_i(s,:)-pp_ris_g,2,2),2*pi))/(norm(pp_scr_i(s,:)-p_bs)*norm(pp_scr_i(s,:)-p_ris))^eta_i;
end

% BS-MU steering vectors:
for kk=1:K

    d_BStoMU = (p_mu(kk,:)-p_bs)/norm(p_mu(kk,:)-p_bs);
    a_bs_mu{kk}(:,1) = exp(1j*kappa*pp_bs*d_BStoMU');

    a_ris_mu{kk}(:,1) = exp(1j*mod(kappa*vecnorm(p_mu(kk,:)-pp_ris_g,2,2),2*pi));

    for s=1:S

         d_BStoScr = (pp_scr_d{kk}(s,:)-p_bs)/norm(pp_scr_d{kk}(s,:)-p_bs);
         a_bs_mu{kk}(:,s+1) = beta*exp(1j*kappa*pp_bs*d_BStoScr')/(norm(pp_scr_d{kk}(s,:)-p_bs)*norm(pp_scr_d{kk}(s,:)-p_mu(kk,:)))^eta_d;

         a_ris_mu{kk}(:,1+s) = beta*exp(1j*mod(kappa*vecnorm(pp_scr_r{kk}(s,:)-pp_ris_g,2,2),2*pi))/(norm(pp_scr_r{kk}(s,:)-p_mu(kk,:))*norm(pp_scr_r{kk}(s,:)-p_ris))^eta_r;
    end
    
end


% Hi LOS channel:
H_i_los = factor_Hi*sqrt(hbar_pl_i)*a_ris_bs(:,1)*a_bs_ris(:,1)';
H_i_nlos = factor_Hi*a_ris_bs(:,2:S)*a_bs_ris(:,2:S)';


% Hr LOS channel:
for kk=1:K
    H_r_los{kk} = factor_Hr*sqrt(hbar_pl_r(kk))*a_ris_mu{kk}(:,1)';
    H_r_nlos{kk} = factor_Hr*sum(a_ris_mu{kk}(:,2:end)');

    H_d_los{kk} = factor_Hd*sqrt(h_blk*hbar_pl_d(kk))*a_bs_mu{kk}(:,1)';
    H_d_nlos{kk} = factor_Hd*sqrt(h_blk)*sum(a_bs_mu{kk}(:,2:end)');
end

%% Combined channel:
H_i = sqrt(Krice_i/(1+Krice_i))*H_i_los+sqrt(1/(1+Krice_i))*H_i_nlos;

for kk=1:K
    H_d{kk} = sqrt(Krice_d/(1+Krice_d))*H_d_los{kk}+sqrt(1/(1+Krice_d))*H_d_nlos{kk};
    H_r{kk} = sqrt(Krice_r/(1+Krice_r))*H_r_los{kk}+sqrt(1/(1+Krice_r))*H_r_nlos{kk};
end

%% Other outputs


Param_output=struct('H_r_los',H_r_los,'H_r_nlos',H_r_nlos,...
    'pp_ris',pp_ris,'pp_ris_g',pp_ris_g,'a_bs_irs',a_bs_ris,'a_irs_bs',a_ris_bs);



