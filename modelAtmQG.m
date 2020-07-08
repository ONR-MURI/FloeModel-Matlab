% This script solves the DIMENSIONAL 2-layer QG with equal layers and a
% rigid lid in a doubly-periodic domain. The spatial dimension is
% discretized by the spectral method while the temporal dimension is
% discretized by Additive Runge–Kutta ARk4. @ Q. Deng 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

% Set simulation parameters
N = 100;        % Number of modes in each direction; 4 times
dt = 1000;        % initial time step size
Nt = 10000;       % number of time steps
qlim = 1.7E1;   % if any q > qlim, simulation stops

% Set physical reference parameters
U0 = 10;        % Horizontal wind
hs = 1e6;       % Horizontal scale, in meters
dz = 4e3;       % Veritical scale, in meters
cf = 1.03e-4;   % Coriolis parameter
t0 = 300;       % reference temperator
be = 1.5e-3;    % Be from Sam's paper
gc = 9.8;       % gravity constant
kb = 1.62e-11;  % Dimensional beta = kb^2 as in beta plan; 1.618624e-11
U = 2;
r = 2e-6;                 % Dimensional Ekman friction coefficient
nu = hs^7*U0*5e-15;     % Coefficient of biharmonic vorticity diffusion
kd = cf^2*t0/(gc * be * dz^2);      % Dimensional deformation wavenumber

c = 0.7;%rand(N);         % c for concentration; can be a constant
Ti = 271.15;         %Sea Ice temperature, in kelvin
To = 273.15;         %Sea ocean temperature, in kelvin
Tc = cf/be;          % theta = Tc \partial_z psi
taub = 1e-6;          % relaxation time scale for buoyancy flux, 10 days
Ts = c*Ti+(1-c)*To-271.15;  % surface temperature
Tsf = fft2(Ts);      % fourier modes

% Put parameters in QG equation into a struct
params = struct('U',U,'kd',kd,'kb',kb, 'r',r, 'nu',nu, 'N',N, 'dt',dt,'Tsf',Tsf,'Tc',Tc, 'taub',taub, 'dz',dz, 'hs',hs);

% Set up hyperviscous PV dissipation
k = [0:N/2 -N/2+1:-1]';  % wavenumbers
L = zeros([N N 2]);
for jj=1:N
    for ii=1:N
        kr = sqrt((k(ii)^2+k(jj)^2)/hs^2); % 2D wave number
        L(ii,jj,:) = -nu*kr^8;      % s=4, this is the hyperviscous term
    end
end
clear kr ii jj

% Initialize
t = 0;
% qp(:,:,2) = randn(params.N);
% qp(:,:,2) = qp(:,:,2)-mean(mean(qp(:,:,2)));
% qp(:,:,1) = qp(:,:,2);
% q = 1e-6*fft2(qp);   % scaling \partial_t q = 1e-10

load('q100.mat','q');
load('qp100.mat','qp');
%load('psi100.mat','psi');

if isempty(dir('dataAtmQG')); mkdir('dataAtmQG'); end

dx = 2*hs/N; xx = -hs:dx:hs; % domain [-hs, hs] = [-1e6, 1e6];
[XX,YY] = meshgrid(xx,xx);

% Diagnostics
countDiag = 10; % was 100, Compute diagnostics every countDiag steps
T = zeros(1,Nt/countDiag);
vb = zeros(1,Nt/countDiag);
utz = zeros(N,Nt/countDiag);
ke = zeros(N/2+1,Nt/countDiag); %kinetic energy
ape = zeros(N/2+1,Nt/countDiag); % average potential energy
figg = figure;
fig3 = figure;

psi = CalcStreamFcn(q,params);
psi1 = real(ifft2(psi));
psi0 = psi1(:,:,2);

meanp = zeros(Nt,1);
medianp = zeros(Nt,1);

%% Main loop, time-marching
ind = 0;
for ii=1:Nt
    %disp(ii);
    
    if mod(ii,countDiag)==0 %&& ii>3000
        if any(isnan(q(:))), break, end
        T(ii/countDiag)=t;
        [KE,APE] = Spectrum(q,params);
        ke(:,ii/countDiag) = KE; ape(:,ii/countDiag) = APE;
        [VB,UTZ] = QG_Diagnostics(q,params);
        vb(ii/countDiag) = VB; utz(:,ii/countDiag) = UTZ;
        if mod(ii, 1000)==0
            display(['iteration i = ', num2str(ii), '; time step dt = ',num2str(dt), ', ene = ',num2str(sum(KE+APE))]);
        end
        
        psi = CalcStreamFcn(q,params);
        psi1 = real(ifft2(psi));
        psi = 0.5* (psi1(:,:,1) + psi1(:,:,2)); % barotrophic mode
        %psi = psi1(:,:,2);
        ind = ind+1;
        
        figure(fig3);
        imagesc(psi);
        %caxis manual
        %caxis([-10 10]);
        saveas(fig3,['./dataAtmQG/' num2str(ind,'psi%03.f') '.jpg'],'jpg');
    end
    
    M = 1./(1-.25*dt*L);
    % First stage ARK4
    k0 = RHS_Spectral(q,params);
    l0 = L.*q;
    % Second stage
    q1 = M.*(q+.5*dt*k0+.25*dt*l0);
    k1 = RHS_Spectral(q1,params);
    l1 = L.*q1;
    % Third stage
    q2 = M.*(q+dt*(13861*k0/62500+6889*k1/62500+8611*l0/62500-1743*l1/31250));
    k2 = RHS_Spectral(q2,params);
    l2 = L.*q2;
    % Fourth stage
    q3 = M.*(q+dt*(-0.04884659515311858*k0-0.1777206523264010*k1+0.8465672474795196*k2...
    +0.1446368660269822*l0-0.2239319076133447*l1+0.4492950415863626*l2));
    k3 = RHS_Spectral(q3,params);
    l3 = L.*q3;
    % Fifth stage
    q4 = M.*(q+dt*(-0.1554168584249155*k0-0.3567050098221991*k1+1.058725879868443*k2...
    +0.3033959883786719*k3+0.09825878328356477*l0-0.5915442428196704*l1...
    +0.8101210538282996*l2+0.2831644057078060*l3));
    k4 = RHS_Spectral(q4,params);
    l4 = L.*q4;
    % Sixth stage
    q5 = M.*(q+dt*(0.2014243506726763*k0+0.008742057842904184*k1+0.1599399570716811*k2...
    +0.4038290605220775*k3+0.2260645738906608*k4+0.1579162951616714*l0...
    +0.1867589405240008*l2+0.6805652953093346*l3-0.2752405309950067*l4));
    k5 = RHS_Spectral(q5,params);
    l5 = L.*q5;
    
    % Successful step, proceed to evaluation
    t = t+dt;
    qp = real(ifft2(q+dt*(0.1579162951616714*(k0+l0)+0.1867589405240008*(k2+l2)+...
    0.6805652953093346*(k3+l3)-0.2752405309950067*(k4+l4)+(k5+l5)/4)));
    q = fft2(qp);
     
    if any(abs(qp(:))>qlim)
        fprintf(['qp = ', num2str(max(abs(qp(:)))),'\n']);
        break
    end
        
end

psi = CalcStreamFcn(q,params);
psi1 = real(ifft2(psi));
psi = psi1;
psi_winds = psi(:,:,2); % the second layer winds enteracting with ices
winds = UpdateWinds(psi_winds,hs,hs,dx);

fig3 = figure;
figure(fig3)
fig3=plotAtmWinds(fig3,winds,0);
