close all; clear all;

%% Set Flags

RIDGING=false; 

FRACTURES=false;

PERIODIC=false;

PACKING = false;

WELDING = false;

CORNERS = false;

COLLISION = true;

AVERAGE = false;

ifPlot = true; %Plot floe figures or not?

ifPlotStress = false;

ifPlotStressStrain = false;

%% Initialize model vars

dt=10; %Time step in sec
% h0 = 0.1; %thickness of ice that gets packed in

%Define ocean currents
nDTpack = 100;
% [ocean, HFo, h0]=initialize_ocean(dt,nDTpack);
[ocean, HFo, h0]=initialize_ocean_inertial(dt,nDTpack);

%Define 10m winds
winds=[0 0];

%Define boundaries
c2_boundary=initialize_boundaries();
Ly = max(c2_boundary(2,:));Lx = max(c2_boundary(1,:));
c2_boundary_poly = polyshape(c2_boundary');
min_floe_size = 4*Lx*Ly/20000;%4*Lx*Ly/25000;

%Initialize Floe state
%Floe=initialize_Floe('FloeShapes.mat');
height.mean = 2;
height.delta = 0;
target_concentration = 1;%[1;1;0];
% [Floe, Nb] = initial_concentration_Nares(c2_boundary,target_concentration,height,50,min_floe_size);
[Floe, Nb] = initial_concentration(c2_boundary,target_concentration,height,250,min_floe_size);
% Floe = Floe(1:200);
%load Floe0; Nb = 0;
if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
end
min_floe_size = (4*Lx*Ly-sum(cat(1,Floe(1:Nb).area)))/25000;
%load('PackedFloesFullDomain.mat','Floe');
%Floe= create_packed_domain();

%%

dhdt = 1;

nDTOut=50; %Output frequency (in number of time steps)

nSnapshots=90; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nPar = 6; %Number of workers for parfor
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(nPar);
else
    delete(poolobj);
    parpool(nPar);
end

target_concentration=1;
tStart = tic; 
doInt.flag = true;
doInt.step = 10;

% specify coarse grid size
LxO= 2*max(ocean.Xo);LyO= 2*max(ocean.Yo);
Nx=3; Ny=fix(Nx*LyO/LxO);
xc = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
yc = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
Xc = (xc(1:end-1)+xc(2:end))/2; Yc = -(yc(1:end-1)+yc(2:end))/2;

%initialize dissolved ice at zero
dissolvedNEW=zeros(Ny,Nx);

%Initiailize Eulearian Data
[eularian_data] = calc_eulerian_stress(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
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
Fx = zeros(Ny, Nx); Fy = zeros(Ny, Nx);
Sig = zeros(Ny, Nx);

%% Calc interactions and plot initial state
Floe=Floe(logical(cat(1,Floe.alive)));
[Floe,dissolvedNEW] = floe_interactions_all(Floe, ocean, winds,c2_boundary, dt,HFo,min_floe_size,Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING); % find interaction points
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
while im_num<nSnapshots

%     ocean.Uocn=ocean.U*cos(ocean.fCoriolis*Time)*ones(size(ocean.Vocn));  
%     ocean.Vocn=-ocean.U*sin(ocean.fCoriolis*Time)*ones(size(ocean.Uocn));
    
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
%         T = (0.05*cos(im_num/dt)+1i*0.05*sin(im_num/dt))/1000;
%         Zhat = (-ocean.r +1i*(ocean.fCoriolis+ocean.Sig))./(ocean.Sig.^2-ocean.fCoriolis.^2-ocean.r^2-2*1i*ocean.r*ocean.Sig)*T/ocean.H;
%         ocean.U = real(ifft(Zhat)); ocean.V = imag(ifft(Zhat));
%         ocean.Uocn = ocean.U*ones(size(ocean.Uocn));
%         ocean.Vocn = ocean.V*ones(size(ocean.Vocn));
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
        Zhat = (-ocean.r +1i*(ocean.fCoriolis+Sig))./(Sig.^2-ocean.fCoriolis.^2-ocean.r^2-2*1i*ocean.r*Sig).*T/ocean.H;
        [~, I2] = max(Zhat);
%         U = real(ifft(Zhat(I2))); V = imag(ifft(Zhat(I2)));
%         nu = fix(log10(abs(U)))+1; nv = fix(log10(abs(V)))+1; n = min([nu nv]);
        %Zhat = (-ocean.r +1i*(ocean.fCoriolis+ocean.Sig))./(ocean.Sig.^2-ocean.fCoriolis.^2-ocean.r^2-2*1i*ocean.r*ocean.Sig)*T/ocean.H;
%         ocean.U = U/(10^n); ocean.V = V/(10^n);
        ocean.Uocn = Amp*cos(Sig(I2)*Time)*ones(size(ocean.Uocn)); ocean.Vocn = Amp*sin(Sig(I2)*Time)*ones(size(ocean.Uocn));
    else
        doInt.flag=false;
    end

    if AVERAGE
%         [eularian_data] = calc_eulerian_stress(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
        [eularian_data] = calc_eulerian_stress(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
        SigXX = SigXX+squeeze(eularian_data.stressxx); SigYX = SigYX+squeeze(eularian_data.stressyx);
        SigXY = SigXY+squeeze(eularian_data.stressxy); SigYY = SigYY+squeeze(eularian_data.stressyy);
        Eux = Eux+squeeze(eularian_data.strainux); Evx = Evx+squeeze(eularian_data.strainvx);
        Euy = Euy+squeeze(eularian_data.strainuy); Evy = Evy+squeeze(eularian_data.strainvy);
        U = U+squeeze(eularian_data.u);V = V+squeeze(eularian_data.v);
        dU = dU+squeeze(eularian_data.du);dV = dV+squeeze(eularian_data.dv);
        Fx = Fx+squeeze(eularian_data.force_x);Fy = Fy+squeeze(eularian_data.force_y);
        Sig = Sig+squeeze(eularian_data.stress);
        [DSig1, DSig2, DSigX, DSigY] = Calc_Stress(eularian_data,dt, c2_boundary);
        DivSigX = DivSigX+DSigX; DivSig1 = DivSig1+DSig1;
        DivSigY = DivSigY+DSigY; DivSig2 = DivSig2+DSig2;
    end
    
    if mod(i_step,nDTOut)==0  %plot the state after a number of timesteps
        
        simp = 0;
        floenew = [];
        for ii = 1:length(Floe)
            floe = Floe(ii);
            if length(Floe(ii).c0) > 30
                floe2 = FloeSimplify(Floe(ii));
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
                simp = simp+1;
            end    
        end
        Floe = [Floe floenew];
        
        if WELDING
            weldrate = 0.1;%Set rate at which floes will meld
            A=cat(1,Floe.area);
            if max(A) > Amax
                Amax = max(A);
            end
            FloeOld = Floe;
            Floe = Weld_Floes(Floe,Nb,weldrate,Amax);
        end

        [eularian_data] = calc_eulerian_stress(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC);
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
        
        if ifPlotStressStrain
            fig2 = figure(fig2);
            if AVERAGE
                SigXX = SigXX/fix(nDTOut); SigYX = SigYX/fix(nDTOut);
                SigXY = SigXY/fix(nDTOut); SigYY = SigYY/fix(nDTOut);
                DivSigX = DivSigX/fix(nDTOut); DivSig1 = DivSig1/fix(nDTOut);
                DivSigY = DivSigY/fix(nDTOut); DivSig2 = DivSig2/fix(nDTOut);
                Eux = Eux/fix(nDTOut); Evx = Evx/fix(nDTOut);
                Euy = Euy/fix(nDTOut); Evy = Evy/fix(nDTOut);
                U = U/fix(nDTOut); V = V/fix(nDTOut);
                dU = dU/fix(nDTOut); dV = dV/fix(nDTOut);
                Fx = Fx/fix(nDTOut); Fy = Fy/fix(nDTOut);
                Sig = Sig/fix(nDTOut);
            else
                SigXX = squeeze(eularian_data.stressxx); SigYX = squeeze(eularian_data.stressyx);
                SigXY = squeeze(eularian_data.stressxy); SigYY = squeeze(eularian_data.stressyy);
                Eux = squeeze(eularian_data.strainux); Evx = squeeze(eularian_data.strainvx);
                Euy = squeeze(eularian_data.strainuy); Evy = squeeze(eularian_data.strainvy);
                U = U+squeeze(eularian_data.u);V = V+squeeze(eularian_data.v);
                dU = dU+squeeze(eularian_data.du);dV = dV+squeeze(eularian_data.dv);
                Fx = Fx+squeeze(eularian_data.force_x);Fy = Fy+squeeze(eularian_data.force_x);
                Sig = Sig+squeeze(eularian_data.stress);
            end
%             if im_num==3
%                 xx = 1; xx(1) =[1 2];
%             end
            SigO = Sig;
%             imagesc(Xc,Yc,abs(Sig)); colorbar;%title('$\sigma_{xx}$','interpreter','latex','fontsize',16); colorbar; caxis([0 1e3])
            subplot(2,4,1); imagesc(Xc,Yc,SigXX); hold on; quiver(Xc,Yc,U*1e5,V*1e6,'k','autoscale','on'); set(gca,'DataAspectRatio',[1 1 1]); title('$\sigma_{xx}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e4 5e4])
            subplot(2,4,2); imagesc(Xc,Yc,SigYX); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'DataAspectRatio',[1 1 1]); title('$\sigma_{yx}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e4 5e4])
            subplot(2,4,5); imagesc(Xc,Yc,SigXY); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'DataAspectRatio',[1 1 1]); title('$\sigma_{xy}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e4 5e4])
            subplot(2,4,6); imagesc(Xc,Yc,SigYY); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'DataAspectRatio',[1 1 1]); title('$\sigma_{yy}$','interpreter','latex','fontsize',16); colorbar; caxis([-5e4 5e4])
            subplot(2,4,3); imagesc(Xc,Yc,abs(Eux)); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'DataAspectRatio',[1 1 1]); title('$E_{11}$','interpreter','latex','fontsize',16); colorbar; caxis([0 1e-8])
            subplot(2,4,4); imagesc(Xc,Yc,abs(Evx)); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'DataAspectRatio',[1 1 1]); title('$E_{21}$','interpreter','latex','fontsize',16); colorbar; caxis([0 1e-8])
            subplot(2,4,7); imagesc(Xc,Yc,abs(Euy)); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'DataAspectRatio',[1 1 1]); title('$E_{12}$','interpreter','latex','fontsize',16); colorbar; caxis([0 1e-8])
            subplot(2,4,8); imagesc(Xc,Yc,abs(Evy)); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'DataAspectRatio',[1 1 1]); title('$E_{22}$','interpreter','latex','fontsize',16); colorbar; caxis([0 1e-8])
            saveas(fig2,['./figs/' num2str(im_num,'Stress%03.f') '.jpg'],'jpg');
            fig3 = figure(fig3);
            subplot(2,2,1); imagesc(Xc,Yc,DivSigX); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'DataAspectRatio',[1 1 1]); title('$\nabla \cdot \sigma_{x} 1$','interpreter','latex','fontsize',16); colorbar; caxis([-0.1 0.1])
            subplot(2,2,2); imagesc(Xc,Yc,DivSigY); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'DataAspectRatio',[1 1 1]); title('$\nabla \cdot \sigma_{y} 1$','interpreter','latex','fontsize',16); colorbar; caxis([-0.1 0.1])
            subplot(2,2,3); imagesc(Xc,Yc,DivSig1); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'DataAspectRatio',[1 1 1]); title('$\nabla \cdot \sigma_{x} 2$','interpreter','latex','fontsize',16); colorbar; caxis([-0.1 0.1])
            subplot(2,2,4); imagesc(Xc,Yc,DivSig2); hold on; quiver(Xc,Yc,U*1e6,V*1e6,'k','autoscale','on'); set(gca,'DataAspectRatio',[1 1 1]); title('$\nabla \cdot \sigma_{y} 2$','interpreter','latex','fontsize',16); colorbar; caxis([-0.1 0.1])
            saveas(fig3,['./figs/' num2str(im_num,'DivStress%03.f') '.jpg'],'jpg');
        end
        
        save(['./Floes/Floe' num2str(im_num,'%07.f') '.mat'],'Floe','eularian_data');
        SigXX = zeros(Ny, Nx); SigYX = zeros(Ny, Nx);
        SigXY = zeros(Ny, Nx); SigYY = zeros(Ny, Nx);
        DivSigX = zeros(Ny, Nx); DivSig1 = zeros(Ny, Nx);
        DivSigY = zeros(Ny, Nx); DivSig2 = zeros(Ny, Nx);
        Eux = zeros(Ny, Nx); Evx = zeros(Ny, Nx);
        Euy = zeros(Ny, Nx); Evy = zeros(Ny, Nx);
        U = zeros(Ny, Nx); V = zeros(Ny, Nx);
        dU = zeros(Ny, Nx); dV = zeros(Ny, Nx);
        Fx = zeros(Ny, Nx); Fy = zeros(Ny, Nx);
        Sig = zeros(Ny, Nx);
        
        M = cat(1,Floe.mass);
        Mtot(im_num) = sum(M)+sum(Vdnew(:));
        
        im_num=im_num+1;  %image number for saving data and coarse vars;
    end

    if PACKING && h0 > 0
        if mod(i_step+1,nDTpack)==0
            height.mean = h0;
            height.delta = 0;
            [Floe,Vd] = pack_ice(Floe,c2_boundary,dhdt,Vd,target_concentration,ocean,height,min_floe_size,PERIODIC);
        end
    end
    
    %Calculate forces and torques and intergrate forward
    [Floe,dissolvedNEW] = floe_interactions_all(Floe, ocean, winds, c2_boundary, dt, HFo,min_floe_size, Nx,Ny,Nb, dissolvedNEW,doInt,COLLISION, PERIODIC, RIDGING);
    
    if FRACTURES && im_num>3 && mod(i_step,10)==0
        overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
        keep=rand(length(Floe),1)>3*overlapArea;
        keep(1:Nb) = ones(Nb,1);
        fracturedFloes=fracture_floe(Floe(~keep),5);
%         fracturedFloes=fracture_leads(Floe(~keep),Nx,Ny,c2_boundary,eularian_data);
        if ~isempty(fracturedFloes)
            Floe=[Floe(keep) fracturedFloes];
        end
%         [Floe] = FracLeads(Floe,Ny,Nx,Nb,c2_boundary,eularian_data);
%         [Floe] = FracIso(Floe,Ny,Nx,Nb,c2_boundary,SigO);
    end
    
    if CORNERS
        stress = zeros(length(Floe),1);
        for ii = 1:length(Floe)
            stress(ii) = trace(abs(Floe(ii).Stress));
        end
        if max(stress)>0
            stresses=stress/max(stress);
        else
            stresses = stress;
        end
        keep=stresses<4*rand(length(Floe),1);
        keep(1:Nb) = ones(Nb,1);
        fracturedFloes=corners(Floe(~keep),Nb);
        if ~isempty(fracturedFloes)
            Floe=[Floe(keep) fracturedFloes];
        end
    end    
    
    %Advect the dissolved mass
    Area=cat(1,Floe.area);
    dissolvedNEW = calc_dissolved_mass(Floe(Area<min_floe_size),Nx,Ny,c2_boundary_poly)+dissolvedNEW;
%     Vdnew = Advect_Dissolved_Ice(Vd,coarseMean,im_num,dissolvedNEW,c2_boundary,dt);
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

end
tEnd = toc(tStart)
%%
