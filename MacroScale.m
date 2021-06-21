clear all, close all

Nx=10; Ny=1; Na = Nx-1;
Lx=6.5e4; Ly=1.5e4;
x=[-1 -1 1 1 -1]*Lx;
y=[-1 1 1 -1 -1]*Ly;
c2_boundary = [x; y];
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
Xc =  (xc(1:end-1)+xc(2:end))/2; Yc = 0;
dx = xc(2)-xc(1); y = [-65e4 65e4];

uNew = zeros(1,Nx+1); v = zeros(1,Nx+1);
uO = uNew+1; SigIce11 = zeros(1,Nx);
SigIce11x = diff(SigIce11)./diff(Xc);
box = (2*65e4)^2/2;
h = 0.5;
height.mean = h; height.delta = 0;
dtmax = 1000;
nPar = 6;

rho_ice=920; % kg/m3

% ice-air drag coefficient
rho_air=1.2;
Cd_atm=1e-3;
U10 = -10; V10 = 0;
Fx_atm=rho_air*Cd_atm*sqrt(U10^2+V10^2)*U10;
tau_atm = Fx_atm/(h*rho_ice);

min_floe_size = 4*Lx*Ly/20000;%4*Lx*Ly/25000;
xx=[-1 -1 1 1 -1]*dx/2; 
yy=[-1 1 1 -1 -1]*Ly;
c2_boundary = [xx; yy];
target_concentration = 1;
[Floe, Nb] = initial_concentration(c2_boundary,target_concentration,height,100,min_floe_size);
for ii = 1:length(xc(2:end-1))
    Floes{ii} = Floe;
    mass(ii) = sum(cat(1,Floe.mass));
    boundaries(ii) = polyshape(c2_boundary');
end
boundariesI = boundaries;
Foa = Fx_atm*area(boundaries);
F = zeros(1,Nx-1); %Fplus = Fminus;
m = mass(1)*ones(1,Nx-1);
SigIce11 = F;
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

istep = 3; time = 0; dt = 1000;
while time<4*dt
    Abound = area(boundaries);
    if istep>10
        xx = 1; xx(1) =[1 2];
    end
%     for ii = 1:length(uNew(2:end-1))
%         if ii<Nx/2
%             accel(ii) = sum(F(1:2*ii))./sum(mass(1:2*ii))+Fx_atm*sum(Abound(1:2*ii))./sum(mass(1:2*ii));
%         elseif ii>Nx/2
%             accel(ii) = sum(F(2*ii+1-Nx:end))./sum(mass(2*ii+1-Nx:end))+Fx_atm*sum(Abound(2*ii+1-Nx:end))./sum(mass(2*ii+1-Nx:end));
%         else
%             accel(ii) = sum(F(1:2*ii))./sum(mass(1:2*ii))+Fx_atm*sum(Abound(1:2*ii))./sum(mass(1:2*ii));
%         end
%     end
%     accel = (F(1:end-1)+F(2:end))./(mass(1:end-1)+mass(2:end))+(Fx_atm*(Abound(1:end-1)+Abound(2:end)))./(mass(1:end-1)+mass(2:end));

    accel = F./m+(Fx_atm*Abound)./(2*m);
    uNew(2:end-1) = uO(2:end-1) + dt*accel;

    %     if max(abs(dtmax*(accel)))>0.025
%         dt = min(abs(0.05./(accel)));
% %     if max(abs(dtmax*(SigIce11x+tau_atm)))>0.05
% %         dt = min(abs(0.05./(SigIce11x+tau_atm)));
%         dt = dt-mod(dt,20);
%         if dt<40
%             dt = 40;
%         end
%         uNew(2:end-1) = uO(2:end-1) + dt*accel;
% %         uNew(2:end-1) = uO(2:end-1) + dt*(SigIce11x+tau_atm);
%     else
%         uNew(2:end-1) = uO(2:end-1) + dtmax*accel;
% %         uNew(2:end-1) = uO(2:end-1) + dtmax*(SigIce11x+tau_atm);
%         dt = dtmax;
%     end
    
    close all; clear fig; clear fig2;
    fig = figure;
    subplot(1,2,1)
    plot(xc,uO,'kx','linewidth',2);  xlabel('X','FontSize',24);  ylabel('$u $','interpreter','latex','fontsize',24); ylim([-0.3 0.3]);
%    fig = imagesc(Xc,Yc,SigIce11x); hold on; quiver(Xc,Yc,u_iceN*1e6,v*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$(\sigma_{11})_x ~ ii$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
%     saveas(fig,['./figs/' num2str(istep,'Macro%03.f') '.jpg'],'jpg');
    subplot(1,2,2)
    plot(xc(2:end-1),SigIce11,'kx','linewidth',2);  xlabel('X','FontSize',24);  ylabel('$(\sigma_{11}) $','interpreter','latex','fontsize',24);  ylim([-2e5 0]); 
%    fig = imagesc(Xc,Yc,SigIce11x); hold on; quiver(Xc,Yc,u_iceN*1e6,v*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$(\sigma_{11})_x ~ ii$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
    saveas(fig,['./figs/' num2str(istep,'MacroSubplots%03.f') '.jpg'],'jpg');
    save(['./FloesSpice/Macro' num2str(istep,'%07.f') '.mat'],'Floes', 'SigIce11x', 'uO','SigIce11','dt');
    
    SigIce = zeros(1,Nx-1);
    m = zeros(1,Nx-1);
    Fx = zeros(1,Nx-1);
%     if istep >1
    
    load(['./FloesSubzero/Floe' num2str(istep,'%07.f') '.mat'],'U','mass','SigXXa');
%     SigXX = (mass(2:2:length(mass)-2).*SigXXa(2:2:length(mass)-2)+mass(3:2:length(mass)-1).*SigXXa(2:2:length(mass)-1))./(mass(2:2:length(mass)-2)+mass(3:2:length(mass)-1));    
    Ua = (mass(1:2:length(mass)-1).*U(1:2:length(mass)-1)+mass(2:2:length(mass)).*U(2:2:length(mass)))./(mass(1:2:length(mass)-1)+mass(2:2:length(mass)));    
%     [Sig,~,Floe,~, bound] = Subzero_Microscale(Ua(5:6),6,dt,Floe,boundary);
%     boundary = polyshape(bound');
%     save(['./FloesSpice/Macro' num2str(ii,'%07.f') '.mat'],'Floe', 'Sig');

    for jj = 1:length(Ua)-1
%         [Sig,Mtot,floe,force, bound] = Subzero_Microscale(uNew(jj:jj+1),6,dt,Floe,boundariesI(jj));
        [Sig,Mtot,floe,force, bound] = Subzero_Microscale(Ua(jj:jj+1),6,dt,Floes{jj},boundaries(jj));
%         [Sig,mass,floe,bound] = Subzero_Microscale(uNew(jj:jj+1),6,dt,Floes{jj},boundaries(jj));
%         if  jj == 2
%         end
        SigIce(jj) = Sig;
        Fx(jj) = force;
        m(jj) = Mtot;
%         DivSig(jj:jj+1) = Fx;
        Floes{jj} = floe;
        boundaries(jj) = polyshape(bound');
    end
    xx = 1; xx(1) =[1 2];
%    m = mass;%(mass(2:2:length(mass)-2)+mass(3:2:length(mass)-1));
%    F = (Fx(1:2:length(Fx)-1)+Fx(2:2:length(Fx)));
    F = Fx;%(Fx(2:2:length(Fx)-2)+Fx(3:2:length(Fx)-1));
    SigIce11 = SigIce;%(mass(2:2:length(mass)-2).*SigIce(2:2:length(mass)-2)+mass(3:2:length(mass)-1).*SigIce(3:2:length(mass)-1))./(SigIce(2:2:length(mass)-2)+mass(3:2:length(mass)-1));
%     area = (areaF(2:2:length(areaF)-2)+areaF(3:2:length(areaF)-1));
%     a = (dU(2:2:length(dU)-2)+dU(3:2:length(dU)-1))/2;
    
    FloesO = Floes;
%     SigIce11 = (Sig11a(2:2:2*(length(Xc)-1)).*mass(2:2:2*(length(Xc)-1))+Sig11a(3:2:2*length(Xc)-1).*mass(3:2:2*length(Xc)-1))./(mass(2:2:2*(length(Xc)-1))+mass(3:2:2*length(Xc)-1));
%     SigIce11x = diff(SigIce11)./diff(Xc);
%     SigIce11x = (DivSig(1:end-1)+DivSig(2:end))/2;
    uO = uNew;
    U(istep,1:length(uNew)) = uNew;
    time = time + dt;
    istep = istep +1;
end

