clear all, close all

Nx=9; Ny=1; Na = Nx-1;
Lx=6.5e4; Ly=1.5e4;
x=[-1 -1 1 1 -1]*Lx;
y=[-1 1 1 -1 -1]*Ly;
c2_boundary = [x; y];
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
Xc =  (xc(1:end-1)+xc(2:end))/2; Yc = 0;
dx = xc(2)-xc(1); y = [-65e4 65e4];

uNew = zeros(1,Nx+1); v = zeros(1,Nx+1);
uO = uNew; SigIce11 = zeros(1,Na);
SigIce11x = diff(SigIce11)./diff(xc(2:end-1));
box = (2*65e4)^2/2;
h = 2;
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
[Floe, Nb] = initial_concentration(c2_boundary,target_concentration,height,75,min_floe_size);
for ii = 1:length(Xc)
    Floes{ii} = Floe;
    m(ii) = sum(cat(1,Floe.mass));
    boundaries(ii) = polyshape(c2_boundary');
end
boundariesI = boundaries;
Foa = Fx_atm*area(boundaries);
Fminus = zeros(Na,1); Fplus = Fminus;
SigIce11 = Fplus;
%% Time Stepping

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(nPar);
else
    delete(poolobj);
    parpool(nPar);
end

istep = 1; time = 0; dt = 100;
while time<24*3600
    A = zeros(Na); D = A; F = zeros(Na,1); Fatm = F;
    A(1,1) = m(1)/(2*m(2)); A(1,2) = 1/2;
    D(1,1)= -1/dx; D(1,2) = 1/dx;
    F(1) = ((m(2) - m(1))*Fplus(1) + (m(1)-m(2))*Fminus(1))/(2*m(1)*m(2));
    Fatm(1) = (Foa(1) + Foa(2))/(2*m(2));
    for ii = 2:Na-1
        A(ii,ii-1) = 1/4; D(ii,ii-1) = -1/(2*dx);
        A(ii,ii) = (m(ii+1)+m(ii-1))*m(ii)/(2*m(ii+1)*m(ii-1));
        A(ii,ii+1) = 1/4; D(ii,ii+1) = 1/(2*dx);
        F(ii) = (m(ii+1)*Fplus(ii) + m(ii-1)*Fminus(ii))/(2*m(ii+1)*m(ii-1));
        Fatm(ii) = (m(ii+1)*Foa(ii-1) + 2*(m(ii+1)+m(ii-1))*Foa(ii)+m(ii-1)*Foa(ii+1))/(4*m(ii+1)*m(ii-1));
    end
    A(Na,Na-1) = 1/2; A(Na,Na) = m(Na)/(2*m(Na-1));
    D(Na,Na-1)= -1/dx; D(Na,Na) = 1/dx;
    F(Na) = ((m(Na) - m(Na-1))*Fplus(Na) + (m(Na-1)-m(Na))*Fminus(Na))/(2*m(Na-1)*m(Na));
    Fatm(Na) = (Foa(Na-1) + Foa(Na))/(2*m(Na-1));
    b = 1/rho_ice*D*SigIce11+F+Fatm;
    accel = A\b;
    uNew(2:end-1) = uO(2:end-1) + dt*accel';
%     xx = 1; xx(1) =[1 2];

%     if max(abs(dtmax*(accel)))>0.05
%         dt = min(abs(0.05./(accel)));
% %     if max(abs(dtmax*(SigIce11x+tau_atm)))>0.05
% %         dt = min(abs(0.05./(SigIce11x+tau_atm)));
%         dt = dt-mod(dt,20);
%         if dt<40
%             dt = 40;
%         end
%         uNew(2:end-1) = uO(2:end-1) + dt*accel';
% %         uNew(2:end-1) = uO(2:end-1) + dt*(SigIce11x+tau_atm);
%     else
%         uNew(2:end-1) = uO(2:end-1) + dtmax*accel';
% %         uNew(2:end-1) = uO(2:end-1) + dtmax*(SigIce11x+tau_atm);
%         dt = dtmax;
%     end
    
    close all; clear fig; clear fig2;
    fig = figure;
    subplot(1,2,1)
    plot(xc(2:end-1),uO(2:end-1),'kx','linewidth',2);  xlabel('X','FontSize',24);  ylabel('$(\sigma_{11})_x $','interpreter','latex','fontsize',24); ylim([-0.35 0]);
%    fig = imagesc(Xc,Yc,SigIce11x); hold on; quiver(Xc,Yc,u_iceN*1e6,v*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$(\sigma_{11})_x ~ ii$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
    subplot(1,2,2)
    plot(xc(2:end-1),SigIce11,'kx','linewidth',2);  xlabel('X','FontSize',24);  ylabel('$(\sigma_{11}) $','interpreter','latex','fontsize',24);  ylim([-1.5e4 0]); 
%    fig = imagesc(Xc,Yc,SigIce11x); hold on; quiver(Xc,Yc,u_iceN*1e6,v*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$(\sigma_{11})_x ~ ii$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
    saveas(fig,['./figs/' num2str(istep,'MacroSubplots%03.f') '.jpg'],'jpg');
    save(['./FloesSpice/Macro' num2str(istep,'%07.f') '.mat'],'Floes', 'SigIce11x', 'uNew','SigIce11','dt');
   
    Sig11a = zeros(1,2*length(Xc));
    mass = Sig11a;
%     if istep >1
        
    for jj = 1:length(Xc)
        [Sig,Mtot,floe,force, bound] = Subzero_Microscale(uNew(jj:jj+1),6,dt,Floe,boundariesI(jj));
%         [Sig,mass,floe,bound] = Subzero_Microscale(uNew(jj:jj+1),6,dt,Floes{jj},boundaries(jj));
%         if Sig > 0 && jj == 9
%             xx = 1; xx(1) =[1 2];
%         end
        Sig11a((2*jj-1):(2*jj)) = Sig;
        ForceX(jj) = sum(force.*Mtot)/sum(Mtot);
        mass((2*jj-1):(2*jj)) = Mtot;
%         DivSig(jj:jj+1) = Fx;
        Floes{jj} = floe;
        boundaries(jj) = polyshape(bound');
    end
    FloesO = Floes;
    Fplus = ForceX(2:Nx);
    Fminus = ForceX(1:Nx-1);
    SigIce11 = (Sig11a(2:2:2*(length(Xc)-1)).*mass(2:2:2*(length(Xc)-1))+Sig11a(3:2:2*length(Xc)-1).*mass(3:2:2*length(Xc)-1))./(mass(2:2:2*(length(Xc)-1))+mass(3:2:2*length(Xc)-1));
    SigIce11x = diff(SigIce11)./diff(xc(2:end-1));
%     SigIce11x = (DivSig(1:end-1)+DivSig(2:end))/2;
    uO = uNew;
    SigIce11 = SigIce11';
    U(istep,1:length(uNew)) = uNew;
    time = time + dt;
    istep = istep +1;
end

