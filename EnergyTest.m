
Nx=6; 
Lx=1e5; Ly=1e5;
x=[-1 -1 1 1 -1]*Lx;
y=[-1 1 1 -1 -1]*Ly;
c2_boundary = [x; y];
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
Xc =  (xc(1:end-1)+xc(2:end))/2; Yc = 0;
dx = xc(2)-xc(1); y = [-65e4 65e4];
rho_ice=920; % kg/m3
dt = 500;

min_floe_size = 4*Lx*Ly/20000;%4*Lx*Ly/25000;
xx=[-1 -1 1 1 -1]*dx; 
yy=[-1 1 1 -1 -1]*Ly;
c2_boundary = [xx; yy];
target_concentration = 1;
height.mean = 2;
height.delta = 0; %max difference between a flow thickness and the mean floe value
FloeNum = 100;%[75 100 125 150 250];
% mass = 0.5*6.6238e+13*ones(1,20);
load FloeEn;
% [FloeEn, Nb] = initial_concentration(c2_boundary,target_concentration,height,FloeNum,min_floe_size);
for ii = 1:length(FloeEn)
    floe2 = FloeSimplify(FloeEn(ii),0);
    if length(floe2) >1
        xx = 1; xx(1) =[1 2];
    end
    FloeEn(ii) = floe2;
    FloeEn(ii).Xi = 0; FloeEn(ii).Yi = 0;
    poly(ii) = polyshape(FloeEn(ii).c0');
end
Kfac = [1000 2000 10000];
tstep = [1 5 10];
rot1 = [0 30 60];
rot2 = [0 45 90];

save('Floe0.mat','FloeEn');

global Modulus;
%load FloeEn;
Modulus = 1.5e3*(mean(sqrt(cat(1,FloeEn.area)))+min(sqrt(cat(1,FloeEn.area))));
save('modulus.mat','Modulus');
clear Kinitial; clear Kfinal
Kinitial = cell(length(Kfac),length(tstep));
Kfinal = Kinitial;
Time = cell(length(Kfac),length(tstep));
%% Loop
% ii = 50;
for jj = 1:length(rot1)
    for kk = 1:length(rot2)
%         jj = 1; kk = 2; ii = 66;
        K0 = 0; KF = 0; p0 = 0; pv = 0; pu = 0; U1 = 0; U2 = 0; V1 = 0; V2 = 0; omeg1 = 0; omeg2 = 0;
        save('Energy5.mat','K0','KF','U1','U2','V1','V2','omeg1','omeg2','pu','pv','p0');
        for ii = 1:length(FloeEn)-1
%             tmod = tstep(kk); kmod = Kfac(jj);
%             En = ii; save('En.mat','En','kmod','tmod');
            r1mod = rot1(kk); r2mod = rot2(jj);
            En = ii; save('En.mat','En','r1mod','r2mod');
            [K,T] = SubzeroET();
%             if (K(end)-K(1))/K(1) >0.05
%                 xx = 1; xx(1)=[1 2];
%             end
        end
        load Energy5;
        Kinitial{jj,kk} = K0;
%         Time{jj,kk} = T;
        Kfinal{jj,kk} = KF;
        save('KE.mat','Kinitial','Kfinal')
    end
end
RS = 1-sum((K0-KF).^2)/sum(K0.^2);

%% 
close all
figure; hold on;
count = 1;
Final = []; Initial =[];
for jj = 1:3
    for kk = 1:3
            KF= Kfinal{jj,kk};
            K0 = Kinitial{jj,kk};
            Initial = [Initial K0];
            Final = [Final KF];
            %semilogx(K0,(KF-K0)./K0,'kx','linewidth',2)
            [I1, I2] = max((KF-K0)./K0);
            if I1 >0.05
                [I2 jj kk]
            end
            %        K0(count) = K(1); KF(count) = K(end);
            count = count +1;
    end
end
        count = 1;
for jj = 1:length(Kfac)
    for kk = 1:length(tstep)
%         tmod = tstep(kk); t = 10*tmod:6000;
        K = Kinitial{jj,kk};
        T = Time{jj,kk};
        plot(T,K)
        K0(count) = K(1); KF(count) = K(end);
        count = count+1;
    end
end
figure
semilogx(K0,(KF-K0)./K0,'ro','linewidth',2)
[~, I1] = max((KF-K0)./K0)
[~, I2] = min((KF-K0)./K0)
