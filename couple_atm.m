function [QG,OU,winds, heat_flux_atm] = couple_atm(Lx,Ly,dXa)
%This function provides coupling between the sea ice and the atmosphere
QG = false;
OU = false;

winds.X=-Lx:dXa:Lx;   %large extent in the x-dir for a channel
winds.Y=-Ly:dXa:Ly; 

[Xatm, Yatm]=meshgrid(winds.X,winds.Y);

winds.U = 10*ones(size(Xatm)); winds.V = 10*ones(size(Yatm));

Tice = -20; Tatm = 2*ones(size(Xatm))+rand(size(Xatm));
heat_flux = 7.4*10^(-4)*(Tice-Tatm)/(72); %cm^2/s
heat_flux_atm = heat_flux/100^2;
end

