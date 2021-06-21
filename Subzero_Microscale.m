function [SigXXa,massA,Floe,SigHist,c2_boundary] = Subzero_Microscale(u,nPar,tstep,Floe,c2_boundary_poly,SigTarget)
close all
uleft = (u(1)-u(2))/2; 
uright = (u(2)-u(1))/2; 
ucell = 0;%(u(1)+u(2))/2;
% uleft = u(1); uright = u(2); ucell = 0;
for i = 1:length(Floe)
    Floe(i).Ui = ucell;
    if isnan(ucell)
        xx = 1; xx(1) =[1 2];
    end
end

%% Set Flags

RIDGING=false; 

FRACTURES=false;

PERIODIC=false;

PACKING = false;

WELDING = false;

CORNERS = false;

COLLISION = true;

AVERAGE = true;

RAFTING = false;

ifPlot = false; %Plot floe figures or not?

ifPlotStress = false;

ifPlotStressStrain = false;

%% Initialize model vars
% dt = tstep/20;
% dt(dt>10)=10; %Time step in sec
dt = 10;
height.mean = 0.5;
height.delta = 0;
% h0 = 0.1; %thickness of ice that gets packed in

%Define ocean currents
nDTpack = 500;
rho_ice=920;
% [ocean, HFo, h0]=initialize_ocean(dt,nDTpack);
[ocean, HFo, h0]=initialize_ocean(dt,nDTpack);
h0 = 0.5;

%Define 10m winds
winds=[-10 0];

%Define boundaries
%c2_boundary=initialize_boundaries();
% Lx=65e4/2; Ly=65e4; 
% x=[-1 -1 1 1 -1]*Lx; 
% y=[-1 1 1 -1 -1]*Ly;
%c2_boundary = [x; y];
%c2_boundary_poly = polyshape(c2_boundary');
c2_boundary = [c2_boundary_poly.Vertices' c2_boundary_poly.Vertices(1,:)'];
x = c2_boundary(1,:); y = c2_boundary(2,:);
Lx = max(c2_boundary(1,:)); Ly = max(c2_boundary(2,:));
% x=c2_boundary(1,:)+[uleft*dt uleft*dt uright*dt uright*dt uleft*dt]*50;
c2_boundary = [x; y];
c2_boundary_poly = polyshape(c2_boundary');
c2_border = scale(c2_boundary_poly,2); c2_border = subtract(c2_border, c2_boundary_poly);

% c2_border = polyshape(2*[-Lx -Lx Lx Lx; -Ly Ly Ly -Ly]'); c2_border = subtract(c2_border, c2_boundary_poly);
floebound = initialize_floe_values(c2_border, height);
min_floe_size = 4*Lx*Ly/20000;%4*Lx*Ly/25000;
% uleft = 0; uright = 0;

%Initialize Floe state
%Floe=initialize_Floe('FloeShapes.mat');
target_concentration = 1; Nb = 0;
%[Floe, Nb] = initial_concentration(c2_boundary,target_concentration,height,500,min_floe_size);
%[Floe, Nb] = initial_concentration(c2_boundary,target_concentration,height,75,min_floe_size);
% Floe = Floe(1:200);
%load Floe0; Nb = 0;
% Floe = Floe0;
if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
end
min_floe_size = (4*Lx*Ly-sum(cat(1,Floe(1:Nb).area)))/25000;
global Modulus
load Modulus;%Modulus = 1.0334e+08;
% Modulus = 1.5e3*(mean(sqrt(cat(1,Floe.area)))+max(sqrt(cat(1,Floe.area))));
%load('PackedFloesFullDomain.mat','Floe');
%Floe= create_packed_domain();

%%

dhdt = 1;

nDTOut=tstep/(dt); %Output frequency (in number of time steps)

nSimp = 50;

nBound = 5;

% poolobj = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(poolobj)
%     parpool(nPar);
% else
%     delete(poolobj);
%     parpool(nPar);
% end

target_concentration=1;
tStart = tic; 
doInt.flag = true;
doInt.step = 10;

% specify coarse grid size
LxO= 2*max(ocean.Xo);LyO= 2*max(ocean.Yo);
Nx=1; Ny=1;%fix(Nx*LyO/LxO);
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
yc = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
Xc = (xc(1:end-1)+xc(2:end))/2; Yc = -(yc(1:end-1)+yc(2:end))/2;

%initialize dissolved ice at zero
dissolvedNEW=zeros(Ny,Nx);

%Initiailize Eulearian Data
[eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
Vd = zeros(Ny,Nx,2);
Vdnew=zeros(Ny, Nx);
SigXX = zeros(Ny, Nx); SigYX = zeros(Ny, Nx);
SigXY = zeros(Ny, Nx); SigYY = zeros(Ny, Nx);
DivSigX = zeros(Ny, Nx); DivSig1 = zeros(Ny, Nx);
DivSigY = zeros(Ny, Nx); DivSig2 = zeros(Ny, Nx);
Eux = zeros(Ny, Nx); Evx = zeros(Ny, Nx);
Euy = zeros(Ny, Nx); Evy = zeros(Ny, Nx);
U = zeros(Ny, Nx); V = zeros(Ny, Nx);
dU = zeros(Ny, Nx); dV = zeros(Ny, Nx);
mass = zeros(Ny, Nx); Fy = zeros(Ny, Nx);
Fx = zeros(Ny, Nx);
Sig = zeros(Ny, Nx);

%% Calc interactions and plot initial state
Floe=Floe(logical(cat(1,Floe.alive)));
[Floe,dissolvedNEW] = floe_interactions_all(Floe, floebound, uright, ucell, ocean, winds,c2_boundary, dt,HFo,min_floe_size,Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING); % find interaction points
A=cat(1,Floe.area);
Amax = max(A);

%% Initialize time and other stuff to zero
if isempty(dir('figs')); disp('Creating folder: figs'); mkdir('figs'); end
if isempty(dir('Floes')); disp('Creating folder: Floes'); mkdir('Floes'); end

if ~exist('Time','var')
    Time=0;
    i_step=0;
    im_num=1;
    fig=0;
    fig2=figure('Position',[100 100 1000 500],'visible','on');
    fig3 = figure;
end


%% Solving for floe trajectories
tic;
while Time<tstep+1

    if mod(i_step,10)==0
        
        disp(' ');
        toc
        disp([num2str(i_step) ' timesteps comleted']); 
        numCollisions = calc_collisionNum(Floe);
        sacked = sum(~cat(1, Floe.alive));
        if sacked>0, disp(['sacked floes: ' num2str(sacked)]); end
        disp(['number of collisions: ' num2str(numCollisions)]);
        disp(' ');
        tic
        doInt.flag=true;
    else
        doInt.flag=false;
    end

    if mod(i_step,nSimp)==0
        floenew = [];
        for ii = 1:length(Floe)
            floe = Floe(ii);
            if length(Floe(ii).c0) > 30
                floe2 = FloeSimplify(Floe(ii),0);
                if isfield(floe2,'poly')
                    floe2=rmfield(floe2,{'poly'});
                end
                for jj = 1:length(floe2)
                    if jj == 1
                        
                        Floe(ii) = floe2(jj);
                    else
                        floenew = [floenew floe2(jj)];
                    end
                end
            end
        end
        Floe = [Floe floenew];
    end
    
    if mod(i_step,nDTOut)==0  %plot the state after a number of timesteps
        

        [eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
%         [eularian_data] = calc_eulerian_data(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
        if ifPlot
            [fig] =plot_basic(fig, Time,Floe,ocean,c2_boundary_poly,Nb);
%             [fig] =plot_Nares(fig, Time,Floe,ocean,c2_boundary_poly,Nb);
            saveas(fig,['./figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
        end
        
        if ifPlotStress
            [fig] =plot_basic_stress(fig, Time,Floe,ocean,c2_boundary_poly,Nb);
            saveas(fig,['./figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
        end
        
        if AVERAGE
            SigXXa = SigXX/fix(i_step); SigYXa = SigYX/fix(i_step);
            SigXYa = SigXY/fix(i_step); SigYYa = SigYY/fix(i_step);
            DivSigXa = DivSigX/fix(i_step); DivSig1a = DivSig1/fix(i_step);
            DivSigYa = DivSigY/fix(i_step); DivSig2a = DivSig2/fix(i_step);
            Eux = Eux/fix(i_step); Evx = Evx/fix(i_step);
            Euy = Euy/fix(i_step); Evy = Evy/fix(i_step);
            U = U/fix(i_step); V = V/fix(i_step);
            dU = dU/fix(i_step); dV = dV/fix(i_step);
            massA = mass/fix(i_step); Fy = Fy/fix(i_step);
            FxA = Fx/fix(i_step);
            Sig = Sig/fix(i_step);

        else
            SigXXa = squeeze(eularian_data.stressxx); SigYXa = squeeze(eularian_data.stressyx);
            SigXYa = squeeze(eularian_data.stressxy); SigYYa = squeeze(eularian_data.stressyy);
            Eux = squeeze(eularian_data.strainux); Evx = squeeze(eularian_data.strainvx);
            Euy = squeeze(eularian_data.strainuy); Evy = squeeze(eularian_data.strainvy);
            U = squeeze(eularian_data.u);V = squeeze(eularian_data.v);
            dU = squeeze(eularian_data.du);dV = squeeze(eularian_data.dv);
            massA = squeeze(eularian_data.Mtot);Fy = squeeze(eularian_data.force_y);
            FxA = squeeze(eularian_data.force_x);
            Sig = squeeze(eularian_data.stress);
        end
        if ifPlotStressStrain
            fig2 = figure(fig2);
            SigO = Sig;
            subplot(2,4,1); imagesc(Xc,Yc,SigXXa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\sigma_{xx}$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
            subplot(2,4,2); imagesc(Xc,Yc,SigYXa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\sigma_{yx}$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
            subplot(2,4,5); imagesc(Xc,Yc,SigXYa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\sigma_{xy}$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
            subplot(2,4,6); imagesc(Xc,Yc,SigYYa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\sigma_{yy}$','interpreter','latex','fontsize',16); colorbar; caxis([-1e6 1e6])
            subplot(2,4,3); imagesc(Xc,Yc,Eux); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$E_{11}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e-6 5e-6])
            subplot(2,4,4); imagesc(Xc,Yc,Evx); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$E_{21}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e-6 5e-6])
            subplot(2,4,7); imagesc(Xc,Yc,Euy); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$E_{12}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e-6 5e-6])
            subplot(2,4,8); imagesc(Xc,Yc,Evy); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$E_{22}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e-6 5e-6])
            saveas(fig2,['./figs/' num2str(im_num,'Stress%03.f') '.jpg'],'jpg');
            fig3 = figure(fig3);
            subplot(2,2,1); imagesc(Xc,Yc,DivSigXa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\nabla \cdot \sigma_{x}~Homogenized$','interpreter','latex','fontsize',16); colorbar; caxis([-0.7 0.7])
            subplot(2,2,2); imagesc(Xc,Yc,DivSigYa); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\nabla \cdot \sigma_{y}~Homogenized$','interpreter','latex','fontsize',16); colorbar; caxis([-0.7 0.7])
            subplot(2,2,3); imagesc(Xc,Yc,DivSig1a); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\nabla \cdot \sigma_{x}~Momentum$','interpreter','latex','fontsize',16); colorbar; caxis([-0.7 0.7])
            subplot(2,2,4); imagesc(Xc,Yc,DivSig2a); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'YDir','normal','DataAspectRatio',[1 1 1]); title('$\nabla \cdot \sigma_{y}~Momentum$','interpreter','latex','fontsize',16); colorbar; caxis([-0.7 0.7])
            saveas(fig3,['./figs/' num2str(im_num,'DivStress%03.f') '.jpg'],'jpg');
        end
        
        save(['./Floes/Floe' num2str(im_num,'%07.f') '.mat'],'Floe','eularian_data','SigXXa','SigXYa', 'SigYYa','DivSigXa','DivSigYa','DivSig1a','DivSig2a');
%         SigXX = zeros(Ny, Nx); SigYX = zeros(Ny, Nx);
%         SigXY = zeros(Ny, Nx); SigYY = zeros(Ny, Nx);
%         DivSigX = zeros(Ny, Nx); DivSig1 = zeros(Ny, Nx);
%         DivSigY = zeros(Ny, Nx); DivSig2 = zeros(Ny, Nx);
%         Eux = zeros(Ny, Nx); Evx = zeros(Ny, Nx);
%         Euy = zeros(Ny, Nx); Evy = zeros(Ny, Nx);
%         U = zeros(Ny, Nx); V = zeros(Ny, Nx);
%         dU = zeros(Ny, Nx); dV = zeros(Ny, Nx);
%         mass = zeros(Ny, Nx); Fy = zeros(Ny, Nx);
%         Fx = zeros(Ny, Nx);
%         Sig = zeros(Ny, Nx);
        
        M = cat(1,Floe.mass);
        Mtot(im_num) = sum(M)+sum(Vdnew(:));
        
        im_num=im_num+1;  %image number for saving data and coarse vars;
    end

    if PACKING && h0 > 0
        if mod(i_step-1,nDTpack)==0
            height.mean = h0;
            height.delta = 0;
%             xx = 1; xx(1) =[1 2];
            [Floe,Vd] = pack_ice(Floe,c2_boundary,dhdt,Vd,target_concentration,ocean,height,min_floe_size,PERIODIC);
        end
    end
    
       
    %Calculate forces and torques and intergrate forward
    [Floe,dissolvedNEW] = floe_interactions_all(Floe, floebound, uright, ucell, ocean, winds, c2_boundary, dt, HFo,min_floe_size, Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING, RAFTING);
    
    if AVERAGE
        %         [eularian_data] = calc_eulerian_stress(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
        [eularian_data] = calc_eulerian_stress2(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
        SigXX = SigXX+squeeze(eularian_data.stressxx); SigYX = SigYX+squeeze(eularian_data.stressyx);
        SigXY = SigXY+squeeze(eularian_data.stressxy); SigYY = SigYY+squeeze(eularian_data.stressyy);
        Eux = Eux+squeeze(eularian_data.strainux); Evx = Evx+squeeze(eularian_data.strainvx);
        Euy = Euy+squeeze(eularian_data.strainuy); Evy = Evy+squeeze(eularian_data.strainvy);
        U = U+squeeze(eularian_data.u);V = V+squeeze(eularian_data.v);
        dU = dU+squeeze(eularian_data.du);dV = dV+squeeze(eularian_data.dv);
        Fx = Fx+squeeze(eularian_data.force_x);Fy = Fy+squeeze(eularian_data.force_y);
        a1 = Floe(1).interactions; if isempty(a1); a1 = [0 0]; end
        a2 =  Floe(2).interactions; if isempty(a2); a2 = [0 0]; end
        SigHist(:,i_step+1) = squeeze(eularian_data.stressxx);
        UHist(:,i_step+1) = squeeze(eularian_data.u);
%         if abs(squeeze(eularian_data.u))>0.008
%             xx = 1; xx(1) =[1 2];
%         end
%         Fxx(:,i_step+1) = zeros(8,1);
%         Fxx(1:length(a1(:,2))+length(a2(:,2)),i_step+1) = [a1(:,2); a2(:,2)];
%         Fxx2(i_step+1) = sum([a1(:,2); a2(:,2)]);
        mass = mass+squeeze(eularian_data.Mtot);
        Sig = Sig+squeeze(eularian_data.stress);
        [DSig1, DSig2, DSigX, DSigY] = Calc_Stress(eularian_data,dt, c2_boundary);
        DivSigX = DivSigX+DSigX; DivSig1 = DivSig1+DSig1;
        DivSigY = DivSigY+DSigY; DivSig2 = DivSig2+DSig2;
    end
    
    if WELDING && mod(i_step,25)==0
	weldrate = 100;%Set rate at which floes will meld
	A=cat(1,Floe.area);
        if max(A) > Amax
           Amax = max(A);
        end
        FloeOld = Floe;
        Floe = Weld_Floes(Floe,Nb,weldrate,Amax);
    end

    if FRACTURES && im_num>1 && mod(i_step,10)==0
        clear Stress, clear A;
        for ii = 1:length(Floe)
            Stress(ii) = max(abs(eig(Floe(ii).Stress)));
        end
        [B,TF] = rmoutliers(Stress);
        keep=rand(length(Floe),1)>Stress'/max(Stress(TF));
        keep(1:Nb) = ones(Nb,1);
        fracturedFloes=fracture_floe(Floe(~keep),3);
        if ~isempty(fracturedFloes)
            Floe=[Floe(keep) fracturedFloes];
        end
    end
    
    if CORNERS && mod(i_step,3)==0
        clear Stress, clear A;
        for ii = 1:length(Floe)
            Stress(ii) = max(abs(eig(Floe(ii).Stress)));
        end
        [B,TF] = rmoutliers(Stress);
        keep=rand(length(Floe),1)>2*Stress'/max(Stress(TF));
        keep(1:Nb) = ones(Nb,1);
        fracturedFloes=corners(Floe(~keep),Nb);
        if ~isempty(fracturedFloes)
            Floe=[Floe(keep) fracturedFloes];
        end

    end   
    
    
    %Advect the dissolved mass
    Area=cat(1,Floe.area);
    dissolvedNEW = calc_dissolved_mass(Floe(Area<min_floe_size),Nx,Ny,c2_boundary_poly)+dissolvedNEW;
    Vdnew = Vd(:,:,1)+dissolvedNEW;
    dissolvedNEW=zeros(Ny,Nx);
    Vd(:,:,2) = Vd(:,:,1);
    Vd(:,:,1) = Vdnew;
    
    Area=cat(1,Floe.area);
    if sum(Area<min_floe_size)>0, display(['num of small floes killed:' num2str(sum(Area<min_floe_size))]); end
    Floe=Floe(Area> min_floe_size);
    live = cat(1,Floe.alive);
    Floe(live == 0) = [];
    
    Time=Time+dt; i_step=i_step+1; %update time index
    if mod(i_step,nBound)==0 && Time < 1001
        x=c2_boundary(1,:)+[uleft*dt uleft*dt uright*dt uright*dt uleft*dt]*nBound;
        c2_boundary = [x; y];
        c2_boundary_poly = polyshape(c2_boundary');
        c2_border = scale(c2_boundary_poly,2); c2_border = subtract(c2_border, c2_boundary_poly);
        floebound = initialize_floe_values(c2_border, height);
    end
    
%     if eularian_data.c < 0.999 && Time >999
%         Time = tstep+1;
%     end
    
end
SigXXa = sum(SigHist)/length(SigHist);%sum(SigHist(i_step-1:i_step))/length(i_step-1:i_step);%SigXX/fix(i_step); SigYXa = SigYX/fix(i_step);
% if max(abs(SigTarget))>3e3
%     xx = 1; xx(1) =[1 2];
% end
            SigXYa = SigXY/fix(i_step); SigYYa = SigYY/fix(i_step);
            DivSigXa = DivSigX/fix(i_step); DivSig1a = DivSig1/fix(i_step);
            DivSigYa = DivSigY/fix(i_step); DivSig2a = DivSig2/fix(i_step);
            Eux = Eux/fix(i_step); Evx = Evx/fix(i_step);
            Euy = Euy/fix(i_step); Evy = Evy/fix(i_step);
            U = U/fix(i_step); V = V/fix(i_step);
            dU = dU/fix(i_step); dV = dV/fix(i_step);
            massA = mass/fix(i_step); Fy = Fy/fix(i_step);
            FxA = Fx/fix(i_step);
            Sig = Sig/fix(i_step);
tEnd = toc(tStart)

end

