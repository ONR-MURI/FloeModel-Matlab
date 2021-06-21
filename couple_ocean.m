function [ocean, c2_boundary, heat_flux,h0]=couple_ocean(transport, Lx, Ly,dXo,dt,nDTOut)

ocean.fCoriolis=1.4e-4; % Coriolis parameter.

ocean.turn_angle=15*pi/180; % turning angle between the stress and surface current due to the Ekman spiral; the angle is positive!

% transport=1e3; % horizontal transport, in m^2/s (controls ocean currents) 
% defining ocean currents

% ocean grid;
%dXo=2000; % in meters

%Xo=(-Lx/2-2*dXo):dXo:(Lx/2+2*dXo); 
Xo=-Lx:dXo:Lx;   %large extent in the x-dir for a channel
Yo=-Ly:dXo:Ly; 

[Xocn, Yocn]=meshgrid(Xo,Yo);

x=[-1 -1 1 1 -1]*Lx/2; 
y=[-1 1 1 -1 -1]*Ly/2;
c2_boundary = [x; y];


%defining ocean streamfunction with some eddies
psi_ocean=transport*sin(2*pi*Xocn/Lx).*sin(2*pi*Yocn/Ly); 


figure; imagesc(psi_ocean);
%title('Ocean Streamfunction');

%calculating ocean velocity field 
Uocn=zeros(size(Xocn)); Vocn=zeros(size(Xocn));
Uocn(2:end,:)=-(psi_ocean(2:end,:)-psi_ocean(1:end-1,:))/dXo; 
Vocn(:,2:end)=(psi_ocean(:,2:end)-psi_ocean(:,1:end-1))/dXo;
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
h0 = real(sqrt(-2*dt*heat_flux*nDTOut));

end