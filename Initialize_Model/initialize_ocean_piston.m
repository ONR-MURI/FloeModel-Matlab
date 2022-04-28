function [ocean, heat_flux, h0]=initialize_ocean_piston(dt,nDTOut)

% defining ocean currents
ocean.fCoriolis=1.4e-4; % Coriolis parameter.

ocean.U = 0.5;

ocean.turn_angle=15*pi/180; % turning angle between the stress and surface current due to the Ekman spiral; the angle is positive!

% ocean grid;
dXo=10000; % in meters

Lx = 4e5; kx = pi/Lx; Ly = 4e5; ky = pi/Ly;
%Lx = 150e4; kx = pi/Lx; Ly = 50e4; ky = pi/Ly;
Xo=-Lx:dXo:Lx; Yo=-Ly:dXo:Ly; 
[Xocn, Yocn]=meshgrid(Xo,Yo);

%defining ocean streamfunction with some eddies
transport=5e3; % horizontal transport, in m^2/s (controls ocean currents) 
psi_ocean=transport/1*(sin(4*kx*Xocn).*sin(4*ky*Yocn));%+cos(ocean.fCoriolis*0)*sin(kx*Xocn).*cos(ky*Yocn)); 

%calculating ocean velocity field 
Uocn=zeros(size(Xocn)); Vocn=zeros(size(Xocn));
%Uocn(2:end,:)=-(psi_ocean(2:end,:)-psi_ocean(1:end-1,:))/dXo; 
%Vocn(:,2:end)=(psi_ocean(:,2:end)-psi_ocean(:,1:end-1))/dXo;

% L = 2*max(Xo); W = 2*max(Yo);
% dx = L/20; dy = W/20;
% x = -L/2:dx:L/2-dx; y = -W/2:dy:W/2-dy;
% k0x=2*pi/L; k0y=2*pi/W;
% [ny,nx] = size(Uocn);
% Uo = fft2(Uocn)/(nx*ny); Vo = fft2(Vocn)/(nx*ny);
% Nx2 = 20; Ny2 = 20;
% Uon = [Uo(:,1:Nx2/2+1) Uo(:,nx-Nx2/2+2:nx)];
% Uon = [Uon(1:Ny2/2+1,:); Uon(ny-Ny2/2+2:ny,:)];
% Von = [Vo(:,1:Nx2/2+1) Vo(:,nx-Nx2/2+2:nx)];
% Von = [Von(1:Ny2/2+1,:); Von(ny-Ny2/2+2:ny,:)];
% [k,l]=meshgrid([0:Nx2/2,-Nx2/2+1:-1]*k0x,[0:Ny2/2,-Ny2/2+1:-1]*k0y);

%adding pure divergence 
%Uocn=Uocn-0.1*Xocn/3e4;
%Vocn=Vocn-0.1*Yocn/3e4;

ocean.Xo=Xo;
ocean.Yo=Yo;
ocean.Xocn = Xocn;
ocean.Yocn = Yocn;
ocean.kx = kx;
ocean.ky = ky;
ocean.Uocn=Uocn;
ocean.Vocn=Vocn;
% ocean.Uon = Uon;
% ocean.Von = Von;
% ocean.k = k; ocean.l = l;

Tice = -20; Tocean = 3*ones(size(Xocn));
heat_flux = 7.4*10^(-4)*(Tice-Tocean)/(72); %cm^2/s
heat_flux = heat_flux/100^2; %m^2/s
h0 = real(sqrt(-2*dt*heat_flux*nDTOut));
h0 = mean(h0(:));
%nDTpack = fix(-h0^2/(2*dt*mean(heat_flux(:)))); %frequency with which packing will be run

end
