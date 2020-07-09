function [QG,OU,winds, heat_flux_atm] = couple_atm(Lx,Ly,dXa)
%This function provides coupling between the sea ice and the atmosphere
QG = false;
OU = false;

Xo=-Lx:dXa:Lx;   %large extent in the x-dir for a channel
Yo=-Ly:dXa:Ly; 

[Xocn, Yocn]=meshgrid(Xo,Yo);

winds.U = 10*ones(size(Xocn)); winds.V = 10*ones(size(Yocn));

Tice = -20; Tatm = 2;
heat_flux = 7.4*10^(-4)*(Tice-Tatm)/(72); %cm^2/s
heat_flux_atm = heat_flux/100^2;
end

