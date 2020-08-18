% function Subzero

close all; clear all;

addpath ./dengwirda-inpoly-ebf47d6/ 

%% Set model flags
RIDGING=true; 

FRACTURES=true;

PERIODIC=true;

SUBFLOES = false;

PACKING = false;

WELDING = false;

ifPlot = true; %Plot floe figures or not?

%% Set the floe domain and couple with Ocean and Atmosphere

%Define coupling paramters
dXo = 4e3; dXa = dXo;
transport=1e4; Lx=2e5; Ly=1e5;
% transport=1e4; Lx=1e6; Ly=1e6;
dt=20; %Time step in sec
rho_ice=920;
nDTOut=100; %Output frequency (in number of time steps)

min_floe_size = 1e6; %Set minimum floe size to keep track of

%Define ocean currents
[ocean, c2_boundary,heat_flux.oc,h0]=couple_ocean(transport, Lx, Ly,dXo,dt,nDTOut,0);
c2_boundary_poly=polyshape(c2_boundary(1,:),c2_boundary(2,:));

%Define Atmosphere
[OU, QG, winds, heat_flux.atm]=couple_atm(Lx,Ly,dXa);

%% Initialize the model

%Set initial mean thickness and variance
height.mean = 2;
height.delta = 0.5; %max difference between a flow thickness and the mean floe value

%Specify target floe concentration
target_concentration=1; % could be a vector

%Generate initiial state
Floe = initialize_floe_field(target_concentration, c2_boundary,ocean,SUBFLOES,height, 100, min_floe_size);%
Nb = 0;

%% Initialize model vars

nSnapshots=1000; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nPar = 4; %Number of workerks for parfor

SackedOB = 0; %initialize number of floes sacked for being out of bounds at zero

dhdt = 1; %Sets probability for ice replenishing open space

%% Define Eulerian grid and coarsening factor

% specify coarse grid size
Lx= 2*max(ocean.Xo);Ly= 2*max(ocean.Yo);
Nx=10; Ny=fix(Nx*Ly/Lx);

%Track floe states
Areas = zeros(2,10,nSnapshots);
Thicknesses = zeros(2,10,nSnapshots);
NumFloes = length(Floe);
Atot = repmat(cat(1,Floe.area),nDTOut,1);
htot = repmat(cat(1,Floe.h),nDTOut,1);

%initialize dissolved ice at zero
dissolvedNEW=zeros(Ny,Nx);
Vd = zeros(Ny,Nx,2);

%%  Calc interactions and plot initial state
[Floe,dissolvedNEW, SackedOB] = floe_interactions_all(Floe, ocean, winds,heat_flux,c2_boundary_poly, dt,dissolvedNEW,SackedOB,Nx,Ny, Nb, RIDGING, PERIODIC,SUBFLOES); % find interaction points
Floe=Floe(logical(cat(1,Floe.alive)));


[eularian_data] = calc_eulerian_data(Floe,Nx,Ny,c2_boundary,PERIODIC);

coarseMean=zeros(9,Ny,Nx,nSnapshots);
coarseSnap=zeros(9,Ny,Nx,nSnapshots);
A=cat(1,Floe.area);
Amax = max(A);
SimpMin = @(A) 9*log10(A);%Function to determine when simplificatoin needs to be done based upon number of vertices

%% Initialize time and other stuff to zero
if isempty(dir('figs')); disp('Creating folder: figs'); mkdir('figs'); end
if isempty(dir('Floes')); disp('Creating folder: Floes'); mkdir('Floes'); end

if ~exist('Time','var')
    Time=0;
    i_step=0;
    im_num=1;
    fig=0;
    fig2=0;
%    EulCoarse=zeros(3, length(cCoarse0(:)),nSnapshots); %allocating memory
end


%% Solving for floe trajectories
tic;
gridArea=area(c2_boundary_poly)/Nx/Ny;
Vdnew=zeros(Ny, Nx);
if ifPlot
    fig3 = figure;
end
count = 1;

%Update the winds
if QG
    [psi12,pv12,qp] = AtmQG(0.7,pvold,i_step,1000,N,dt*nDTOut);
    psi1 = real(ifft2(psi12));
    psi_winds = psi1(:,:,2); % the second layer winds enteracting with ices
    winds = UpdateWinds(psi_winds,Lx,Ly,dXa);
    pvold = pv12;
end

while im_num<nSnapshots
     
    display(i_step);
    if mod(i_step,10)==0
        disp(newline);
        elapsedTime = toc;
        disp(num2str(elapsedTime));
        StepTimes(count) = elapsedTime/length(Floe);
        disp([num2str(i_step) ' timesteps comleted']); 
        numCollisions = calc_collisionNum(Floe);
        sacked = sum(~cat(1, Floe.alive));
        if sacked>0, disp(['sacked floes: ' num2str(sacked)]); end
        if SackedOB>0, disp(['total sacked floes for being out of bounds: ' num2str(SackedOB)]); end
        disp(['number of collisions: ' num2str(numCollisions)  newline]);
        tic
        count = count+1;
    end

    %Plot, calculate mean values, and pack new ice after a number of
    %timesteps
    if mod(i_step,nDTOut)==0 
        
        %Update the winds
        if OU
            disp('Running atmosphere OU equations...');
            [tpsi] = OU2D(theta, omega, sigma, mu, psivar, psi, N, dt*nDTOut);
            psi = real(ifft2(tpsi));
            winds = UpdateWinds(psi,Lx,Ly,dXa);
            psi = tpsi;
        elseif QG
            disp('Running atmosphere QG equations...');
            [psi12,pv12,qp] = AtmQG(0.7,pvold,i_step,10,N,dt*nDTOut);
            psi1 = real(ifft2(psi12));
            psi_winds = 0.5* ( psi1(:,:,1) + psi1(:,:,2) ); % barotropic mode
            winds = UpdateWinds(psi_winds,Lx,Ly,dXa);
            pvold = pv12;
        end
            
        ocean.Uocn=U*cos(ocean.fCoriolis*Time);  
        ocean.Vocn=-U*sin(ocean.fCoriolis*Time);
        
        %Floe thickness and area statistics
        [A1,A2] = hist(Atot,10);
        [h1,h2] = hist(htot,10);
        FloeStats(im_num).Num = fix(NumFloes);
        FloeStats(im_num).DissolvedMass = sum(sum(Vdnew));
        FloeStats(im_num).A = Atot;
        FloeStats(im_num).h = htot;
        FloeStats(im_num).Floes = Floe;
        Atot = [];
        htot = [];
        NumFloes = 0;
        save('FloeStats.mat','FloeStats','Amax')
        save(['./Floes/Floe' num2str(im_num,'%07.f') '.mat'],'Floe');
        
        %Run packing function
        if mod(i_step,nDTOut)==0 && PACKING && max(max(h0))>0
            height.mean = h0;
            height.delta = 0;
            [Floe,Vd] = pack_ice(Floe,c2_boundary,dhdt,Vd,target_concentration,ocean,height, min_floe_size, SUBFLOES, PERIODIC);

        end
                
        %Corase mean data
        [eularian_data] = calc_eulerian_data(Floe,Nx,Ny,c2_boundary,PERIODIC);
        coarseSnap(1,:,:,im_num)=eularian_data.c;
        coarseSnap(2,:,:,im_num)=eularian_data.u;
        coarseSnap(3,:,:,im_num)=eularian_data.v;
        coarseSnap(4,:,:,im_num)=eularian_data.du;
        coarseSnap(5,:,:,im_num)=eularian_data.dv;
        coarseSnap(6,:,:,im_num)=eularian_data.mom_x;
        coarseSnap(7,:,:,im_num)=eularian_data.mom_y;
        coarseSnap(8,:,:,im_num)=eularian_data.force_x;
        coarseSnap(9,:,:,im_num)=eularian_data.force_y;        
        save('coarseData.mat','coarseSnap','coarseMean');
        
%         if dhdt > 0
%             c = eularian_data.c;
%             c(c>0.99) = 1;
%             dissolvedNEW = dissolvedNEW + (1-c)*gridArea*rho_ice*mean(mean(h0)); %saying here that open water is being populated by sea ice growth consistent with 0.2 m thick ice
%         end
        
        %Check to see if any floes need to be simplified
        floenew = [];
        for ii = 1:length(Floe)
            ddx = 100;
            if length(Floe(ii).poly.Vertices) > SimpMin(Floe(ii).area)
                floe = FloeSimplify(Floe(ii), ddx,SUBFLOES);
                for jj = 1:length(floe)
                    if jj == 1
                        Floe(ii) = floe(jj);
                    else
                        floenew = [floenew floe(jj)];
                    end
                end
            end
        end
        Floe = [Floe floenew];
        A = cat(1,Floe.area);
        Floe(A<min_floe_size) = [];
        live = cat(1,Floe.alive);
        Floe(live==0)=[];

        %Plot the floes
        if ifPlot
            [fig]=plot_basic(fig, Time,Floe, ocean, c2_boundary_poly);
            %[fig, fig2]=plot_Floes(fig,fig2, Time,Floe, ocean, c2_boundary_poly, im_num,PERIODIC);
            if im_num>1
            figure(fig3);
            imagesc(Vdnew/(gridArea*rho_ice)); axis xy
            u=squeeze(coarseMean(2,:,:,im_num));
            v=squeeze(coarseMean(3,:,:,im_num));
            colormap('gray');colorbar;
            hold on; quiver(u,v,'r')
            drawnow
            set(0,'CurrentFigure',fig);
            end

        end
        
        
        im_num=im_num+1;  %image number for saving data and coarse vars;
    end
    
    %Calculate forces and torques and intergrate forward
    [Floe,dissolvedNEW, SackedOB] = floe_interactions_all(Floe, ocean, winds,heat_flux, c2_boundary_poly, dt,dissolvedNEW,SackedOB,Nx,Ny, Nb, RIDGING, PERIODIC,SUBFLOES);
    
    %Run welding every so many time steps
    if mod(i_step-1,nDTOut)==0
        if WELDING
            weldrate = 0.1;%Set rate at which floes will meld
            A=cat(1,Floe.area);
            if max(A) > Amax
                Amax = max(A);
            end
            Floe = Weld_Floes(Floe,weldrate,Amax,SUBFLOES);
        end

    end
    
    
    %Run fracture fucture
    if FRACTURES
        overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
        keep=rand(length(Floe),1)>overlapArea;
        if Nb > 0
            keep(1:Nb) = ones(Nb,1);
        end
        fracturedFloes=fracture_floe(Floe(~keep),3);
        if ~isempty(fracturedFloes), fracturedFloes=rmfield(fracturedFloes,'potentialInteractions');
            Floe=[Floe(keep) fracturedFloes];
        end
    end
    
    %Calculate coarse data to track mean values
    [eularian_data] = calc_eulerian_data(Floe,Nx,Ny,c2_boundary,PERIODIC);    
    coarseMean(1,:,:,im_num)=squeeze(coarseMean(1,:,:,im_num))+eularian_data.c/nDTOut;
    coarseMean(2,:,:,im_num)=squeeze(coarseMean(2,:,:,im_num))+eularian_data.u/nDTOut;
    coarseMean(3,:,:,im_num)=squeeze(coarseMean(3,:,:,im_num))+eularian_data.v/nDTOut;
    coarseMean(4,:,:,im_num)=squeeze(coarseMean(4,:,:,im_num))+eularian_data.du/nDTOut;
    coarseMean(5,:,:,im_num)=squeeze(coarseMean(5,:,:,im_num))+eularian_data.dv/nDTOut;
    coarseMean(6,:,:,im_num)=squeeze(coarseMean(6,:,:,im_num))+eularian_data.mom_x/nDTOut;
    coarseMean(7,:,:,im_num)=squeeze(coarseMean(7,:,:,im_num))+eularian_data.mom_y/nDTOut;
    coarseMean(8,:,:,im_num)=squeeze(coarseMean(8,:,:,im_num))+eularian_data.force_x/nDTOut;
    coarseMean(9,:,:,im_num)=squeeze(coarseMean(9,:,:,im_num))+eularian_data.force_y/nDTOut;
    
    %Find thickness and height data to keep track of these statistics
    NumFloes = NumFloes+length(Floe)/nDTOut;
    Asnap = cat(1,Floe.area);
    Atot = [Atot; Asnap(Nb+1:end)];
    hsnap = cat(1,Floe.h);
    htot = [htot; hsnap(Nb+1:end)];
    
    %Advect the dissolved mass
    Area=cat(1,Floe.area);
    dissolvedNEW = calc_dissolved_mass(Floe(Area<min_floe_size),Nx,Ny,c2_boundary_poly)+dissolvedNEW;
    Vdnew = Advect_Dissolved_Ice(Vd,coarseMean,im_num,dissolvedNEW,c2_boundary,dt);
    dissolvedNEW=zeros(Ny,Nx);
    Vd(:,:,2) = Vd(:,:,1);
    Vd(:,:,1) = Vdnew;
    
    
    if sum(Area<min_floe_size)>0, display(['num of small floes killed:' num2str(sum(Area<min_floe_size))]); end
    Floe=Floe(Area> min_floe_size);
    Time=Time+dt; i_step=i_step+1; %update time index

end
%%
% end