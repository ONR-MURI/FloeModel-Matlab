clear all, close all

Nx=10; 
Lx=6.5e4; Ly=1.5e4;
x=[-1 -1 1 1 -1]*Lx;
y=[-1 1 1 -1 -1]*Ly;
c2_boundary = [x; y];
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
Xc =  (xc(1:end-1)+xc(2:end))/2; Yc = 0;
dx = xc(2)-xc(1); y = [-65e4 65e4];
rho_ice=920; % kg/m3
dt = 1000;

min_floe_size = 4*Lx*Ly/20000;%4*Lx*Ly/25000;
xx=[-1 -1 1 1 -1]*dx/2; 
yy=[-1 1 1 -1 -1]*Ly;
c2_boundary = [xx; yy];
target_concentration = 1;
height.mean = 2;
height.delta = 0; %max difference between a flow thickness and the mean floe value
[FloeI, ~] = initial_concentration(c2_boundary,target_concentration,height,100,min_floe_size);
boundaryI = polyshape(c2_boundary');
for ii = 1:length(xc(2:end-1))
    Floes{ii} = FloeI;
    boundaries(ii) = polyshape(c2_boundary');
end

nPar = 6;

% xx = 1; xx(1) =[1 2];
%% Time Stepping
% 
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(nPar);
else
    delete(poolobj);
    parpool(nPar);
end

figure; hold on
dt = 1000; 
for istep = 3:10
    SigIce = zeros(1,Nx-1);
    
    load(['./FloesSubzero/Floe' num2str(istep,'%07.f') '.mat'],'U','mass','SigXXa');
    SigXX = (mass(2:2:length(mass)-2).*SigXXa(2:2:length(mass)-2)+mass(3:2:length(mass)-1).*SigXXa(2:2:length(mass)-1))./(mass(2:2:length(mass)-2)+mass(3:2:length(mass)-1));    
    Ua = (mass(1:2:length(mass)-1).*U(1:2:length(mass)-1)+mass(2:2:length(mass)).*U(2:2:length(mass)))./(mass(1:2:length(mass)-1)+mass(2:2:length(mass)));    

%     jj = 5;
    for jj = 1:length(Ua)-1
%         [Sig,Mtot,floe,force, bound] = Subzero_Microscale(uNew(jj:jj+1),6,dt,Floe,boundariesI(jj));
        [Sig,~,floe,~, bound] = Subzero_Microscale(Ua(jj:jj+1),6,dt,Floes{jj},boundaries(jj));
        SigIce(jj) = Sig;
        Floes{jj} = floe;
        boundaries(jj) = polyshape(bound');
    end
%     plot(t,SigXX(5),'kx','linewidth',2);
%     plot(t,Sig,'rx','linewidth',2);
    save(['./FloesSpice/Macro' num2str(istep,'%07.f') '.mat'],'Floes', 'SigIce');    
    FloesO = Floes;
end

close all
t = 3000:1000:75000;
figure, hold on
for ii = 3:75
    load(['./FloesSubzero/Floe' num2str(ii,'%07.f') '.mat']);
    Ua = (mass(1:2:length(mass)-1).*U(1:2:length(mass)-1)+mass(2:2:length(mass)).*U(2:2:length(mass)))./(mass(1:2:length(mass)-1)+mass(2:2:length(mass)));
    Uspice(ii-2,:) = Ua;
    load(['./FloesSubzero/Floe' num2str(ii+1,'%07.f') '.mat']);
    load(['./FloesSpice/Macro' num2str(ii,'%07.f') '.mat']);
    m = (mass(2:2:length(mass)-2)+mass(3:2:length(mass)-1));
    SigXX = (mass(2:2:length(mass)-2).*SigXXa(2:2:length(mass)-2)+mass(3:2:length(mass)-1).*SigXXa(2:2:length(mass)-1))./(mass(2:2:length(mass)-2)+mass(3:2:length(mass)-1));    
    SigSpice(ii-2,:) = SigIce-Sig;
    SigSub(ii-2,:) = SigXX-Sig;
    %     plot(t,SigXX(5),'kx','linewidth',2);
%     plot(t,SigIce(5),'ro','linewidth',2);
%     t = t+dt;
end
