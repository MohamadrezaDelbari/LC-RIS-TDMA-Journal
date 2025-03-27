function [pp_bs,pp_ris] = func_bs_irs_ant_p(Param)

N_bs_x = Param.N_bs_x; % number of BS's antenna along x axis
N_bs_z = Param.N_bs_z; % number of BS's antenna along z axis
N_bs = N_bs_x*N_bs_z;

d_bs_x = Param.d_bs_x; % BS element spacing along x axis
d_bs_z = Param.d_bs_z; % BS element spacing along z axis

N_y = Param.N_y; % number of RIS elements along y axis
N_z = Param.N_z; % number of RIS elements along z axis
N = N_y*N_z;

d_ris_y = Param.d_ris_y; % RIS element spacing along y axis
d_ris_z = Param.d_ris_z; % RIS element spacing along z axis

%% Antenna positions:
% BS antennas
pp_bs_ant_x = linspace(-(N_bs_x-1)*d_bs_x/2,(N_bs_x-1)*d_bs_x/2,N_bs_x);
pp_bs_ant_z = linspace(-(N_bs_z-1)*d_bs_z/2,(N_bs_z-1)*d_bs_z/2,N_bs_z);
[XX,ZZ] = meshgrid(pp_bs_ant_x,pp_bs_ant_z);
vx=reshape(XX,[],1);
vz=reshape(ZZ,[],1);
pp_bs = [vx,zeros(N_bs,1),vz];

% RIS elemnts

pp_ris_y = linspace(-(N_y-1)*d_ris_y/2,(N_y-1)*d_ris_y/2,N_y);
pp_ris_z = linspace(-(N_z-1)*d_ris_z/2,(N_z-1)*d_ris_z/2,N_z);
[YY,ZZ] = meshgrid(pp_ris_y,pp_ris_z);
vy=reshape(YY,[],1);
vz=reshape(ZZ,[],1);
pp_ris = [zeros(N,1),vy,vz];

end