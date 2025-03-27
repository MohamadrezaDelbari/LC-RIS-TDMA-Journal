%%----------------------------------------------------%%
%%----- Mohamadreza Delbari
%%----- Please cite: Fast Reconfiguration of LC-RISs in TDMA Protocol
%%----- Please cite: DOI: 10.1109/ICCWorkshops59551.2024.10615422
%%----------------------------------------------------%%

function [W] = NF_RIS_design(f,Ny,Nz,pr,pt,ppirs)
c=3*10^8;
lamda=c/f;
k=2*pi/lamda;
N=Ny*Nz;
pirs=mean(ppirs);
%______________________Position________________________
prx=linspace(pr(1,1),pr(2,1),Nz); %Receiver
pry=linspace(pr(1,2),pr(2,2),Ny);
% Create a grid of these points
[gridY, gridX] = meshgrid(pry, prx);

% Flatten the grid to vectors
vectorX = gridX(:);
vectorY = gridY(:);

% Concatenate these vectors into a Ny*Nz x 2 matrix
pr_vector = [vectorX, vectorY pr(1,3)*ones(N,1)];

dt=vecnorm([pt(1,1) pt(1,2) pt(1,3)]'-ppirs');
dr=vecnorm(pr_vector'-ppirs');
dn=vecnorm(pirs'*ones(1,N)-[vectorX vectorY pr(1,3)*ones(N,1)]');
%dn=vecnorm(pirs'-[prx' pry' pr(1,3)*ones(N,1)]');

w = exp(1i*k*(dr-dn-dt));
W=diag(w);
%save("near_ana_1D.mat","w","W")
end
