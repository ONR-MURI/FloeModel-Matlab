function [ocean, heat_flux, nDTpack]=initialize_ocean_Nares(dt,h0)

% defining ocean currents
ocean.fCoriolis=1.4e-4; % Coriolis parameter.

ocean.U = 0.5;

ocean.turn_angle=15*pi/180; % turning angle between the stress and surface current due to the Ekman spiral; the angle is positive!

% ocean grid;
dXo=2000; % in meters

dXo = 4e3; dXa = dXo;
Lx=2e5; Ly=1e6;
Xo=-1.1*Lx/2:dXo:1.1*Lx/2; Yo=-1.1*Ly/2:dXo:1.1*Ly/2; 
[Xocn, Yocn]=meshgrid(Xo,Yo);

%defining ocean streamfunction with some eddies
transport=1e3; % horizontal transport, in m^2/s (controls ocean currents) 
psi_ocean=transport*sin(2*pi*Xocn/40e3).*cos(2*pi*Yocn/50e3); 

%calculating ocean velocity field 
Uocn=zeros(size(Xocn)); Vocn=-0.5*ones(size(Xocn));
% Uocn(2:end,:)=-(psi_ocean(2:end,:)-psi_ocean(1:end-1,:))/dXo; 
% Vocn(:,2:end)=(psi_ocean(:,2:end)-psi_ocean(:,1:end-1))/dXo;

%adding pure divergence 
%Uocn=Uocn-0.1*Xocn/3e4;
%Vocn=Vocn-0.1*Yocn/3e4;

ocean.Xo=Xo;
ocean.Yo=Yo;
ocean.Uocn=Uocn;
ocean.Vocn=Vocn;

Tice = -20; Tocean = 2*ones(size(Xocn));
heat_flux = 7.4*10^(-4)*(Tice-Tocean)/(72); %cm^2/s
heat_flux = heat_flux/100^2; %m^2/s
% h0 = real(sqrt(-2*dt*heat_flux*nDTOut));
nDTpack = fix(-h0^2/(2*dt*mean(heat_flux(:)))); %frequency with which packing will be run

end