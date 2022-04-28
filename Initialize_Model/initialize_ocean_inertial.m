function [ocean, heat_flux, h0]=initialize_ocean_inertial(dt,nDTOut)

% defining ocean currents
ocean.fCoriolis=1.4e-4; % Coriolis parameter.
ocean.r = 1/(2.5*24*3600); % Damping parameter
ocean.H = 15; %ocean mix layer depth

ocean.U = 0.5; ocean.V = ocean.U;

ocean.turn_angle=15*pi/180; % turning angle between the stress and surface current due to the Ekman spiral; the angle is positive!

% ocean grid;
dXo=500; % in meters

Xo=-3e4:dXo:3e4; Yo=-3e4:dXo:3e4; 
[Xocn, Yocn]=meshgrid(Xo,Yo);
ocean.Xocn = Xocn; ocean.Yocn = Yocn;
n = 0:0.01:20; % cycles per day
ocean.Sig = ocean.fCoriolis;%n/(24*3600);

% Generate Red Noise
% signal parameters
Amp = 0.5;
fs = 44100;         % sampling frequency, Hz
T = 60;             % signal duration, s
N = round(fs*T);    % number of samples
% noise generation
x = rednoise(1, N); % red noise, PSD falls off by 10 dB/dec
% calculate the noise PSD
winlen = 2*fs;
window = hann(winlen, 'periodic');
noverlap = winlen/2;
nfft = winlen;
[Pxx, f] = pwelch(x, window, noverlap, nfft, fs, 'onesided');
Sig = f/1e5;
Pxx = Pxx(Sig>1e-4)/1e5; Pxx = Pxx(1:100); Sig = Sig(1:100);
T = (Pxx+1i*Pxx)/1000;
%Zhat = (-ocean.r +1i*(ocean.fCoriolis+ocean.Sig))./(ocean.Sig.^2-ocean.fCoriolis.^2-ocean.r^2-2*1i*ocean.r*ocean.Sig)*T/ocean.H;
Zhat = (-ocean.r +1i*(ocean.fCoriolis+Sig))./(Sig.^2-ocean.fCoriolis.^2-ocean.r^2-2*1i*ocean.r*Sig).*T/ocean.H;
[~, I2] = max(Zhat);
% U = real(ifft(Zhat(I2))); V = imag(ifft(Zhat(I2)));
% nu = fix(log10(abs(U)))+1; nv = fix(log10(abs(V)))+1; n = min([nu nv]);
% ocean.U = U/(10^n); ocean.V = V/(10^n);
% ocean.U = Amp*cos(Sig(I2)*0); ocean.V = Amp*sin(Sig(I2)*0);

%defining ocean streamfunction with some eddies
transport=5e3; % horizontal transport, in m^2/s (controls ocean currents) 
psi_ocean=transport*sin(2*pi*Xocn/40e4).*cos(2*pi*Yocn/50e4); 

%calculating ocean velocity field 
Uocn=zeros(size(Xocn)); Vocn=zeros(size(Xocn));
% Uocn = ocean.U*ones(size(Xocn)); 
% Vocn = ocean.V*ones(size(Xocn));


L = 2*max(Xo); W = 2*max(Yo);
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
ocean.Uocn=Uocn;
ocean.Vocn=Vocn;
ocean.Uocn_p=Uocn;
ocean.Vocn_p=Vocn;
% ocean.Uon = Uon;
% ocean.Von = Von;
% ocean.k = k; ocean.l = l;

Tice = -20; Tocean = 2*ones(size(Xocn));
heat_flux = 7.4*10^(-4)*(Tice-Tocean)/(72); %cm^2/s
heat_flux = heat_flux/100^2; %m^2/s
h0 = real(sqrt(-2*dt*heat_flux*nDTOut));
h0 = mean(h0(:));
%nDTpack = fix(-h0^2/(2*dt*mean(heat_flux(:)))); %frequency with which packing will be run

ocean.TauAtmX_p=zeros(size(Xocn));
ocean.TauIceX_p=zeros(size(Xocn));
ocean.TauAtmY_p=zeros(size(Xocn));
ocean.TauIceY_p=zeros(size(Xocn));


end