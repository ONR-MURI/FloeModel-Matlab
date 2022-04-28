function [ocean, heat_flux, h0]=initialize_ocean_Nares(dt,nDTOut)

% defining ocean currents
ocean.fCoriolis=1.4e-4; % Coriolis parameter.

ocean.U = 0;

ocean.turn_angle=15*pi/180; % turning angle between the stress and surface current due to the Ekman spiral; the angle is positive!

% ocean grid;
dXo=2000; % in meters

dXo=20000; % in meters
Lx=2e5; Ly=1e6;
Xo=-1.5*Lx/2:dXo:1.5*Lx/2; Yo=-1.5*Ly/2:dXo:1.5*Ly/2; 
[Xocn, Yocn]=meshgrid(Xo,Yo);
ocean.Xocn = Xocn; ocean.Yocn = Yocn;

%defining ocean streamfunction with some eddies
transport=5e3; % horizontal transport, in m^2/s (controls ocean currents) 
psi_ocean=transport*sin(2*pi*Xocn/40e3).*cos(2*pi*Yocn/50e3); 

%calculating ocean velocity field 
Uocn=zeros(size(Xocn)); Vocn=-0*ones(size(Xocn));%zeros(size(Xocn));%
% Uocn(2:end,:)=-(psi_ocean(2:end,:)-psi_ocean(1:end-1,:))/dXo; 
% Vocn(:,2:end)=(psi_ocean(:,2:end)-psi_ocean(:,1:end-1))/dXo;

%adding pure divergence 
%Uocn=Uocn-0.1*Xocn/3e4;
%Vocn=Vocn-0.1*Yocn/3e4;

ocean.Xo=Xo;
ocean.Yo=Yo;
ocean.Uocn=Uocn;
ocean.Vocn=Vocn;
ocean.Uocn_p=Uocn;
ocean.Vocn_p=Vocn;

Tice = -20; Tocean = 2*ones(size(Xocn));
heat_flux = 7.4*10^(-4)*(Tice-Tocean)/(72); %cm^2/s
heat_flux = heat_flux/100^2; %m^2/s
h0 = real(sqrt(-2*dt*heat_flux*nDTOut));
h0 = mean(h0(:));

ocean.TauAtmX_p=zeros(size(Xocn));
ocean.TauIceX_p=zeros(size(Xocn));
ocean.TauAtmY_p=zeros(size(Xocn));
ocean.TauIceY_p=zeros(size(Xocn));
end