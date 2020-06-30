%% Atmosphere-sea-ice onw-way coupling model: DEM for ice; Two-layer QG for atmosphere; no ocean.

close all; clear all; clc,

%addpath ~/siam0627/dengwirda-inpoly-ebf47d6/
addpath ~/Documents/Projects/2020/siam0627/dengwirda-inpoly-ebf47d6/

%% Initialize model vars
RIDGING=false;
FRACTURES=false;
PERIODIC=true;
SUBFLOES = false;
PACKING = false;

N = 100;
dXo = 2e6/N; dXa = dXo;
transport=1e5; Lx=1e6; Ly=1e6;

%Define ocean currents
[ocean,c2_boundary]=initialize_ocean_Gyre(transport, Lx, Ly,dXo); %(transport, Lx, Ly,dXo)
c2_boundary_poly=polyshape(c2_boundary(1,:),c2_boundary(2,:));
ocean.Uocn = 0*ocean.Uocn; ocean.Vocn =  0*ocean.Vocn;

load('WindsInit100.mat','winds');
load('PVInit100.mat','pvold');
load('FloeInit50.mat','Floe');

%Initialize Floe state
height.mean = 2;
height.delta = 0.5; %max difference between a flow thickness and the mean floe value
dt=100; %was 20 Time step in sec
nDTOut=10; %Output frequency (in number of time steps)
nSnapshots=100; %Total number of model snapshots to save
nDT=nDTOut*nSnapshots; %Total number of time steps
nPar = 4; %Number of workerks for parfor
ifPlot = true; %Plot floe figures or not?
SackedOB = 0; %initialize number of floes sacked for being out of bounds at zero
dhdt = 0.1; %Rate at which ice thickness increases thermodynamically per day
heat_flux = 0.05/(24*3600); %Rate at which ice thickness increases thermodynamically per day

% specify coarse grid size
Nx=10; Ny=fix(Nx*Ly/Lx);
%initialize dissolved ice at zero
dissolvedNEW=zeros(Ny,Nx);

% Calc interactions and plot initial state
[Floe,dissolvedNEW, SackedOB] = floe_interactions_all_doublePeriodicBCs_bpm(Floe,ocean,winds,heat_flux,c2_boundary_poly,dt,dissolvedNEW,SackedOB,Nx,Ny,RIDGING,PERIODIC,SUBFLOES); % find interaction points
Floe=Floe(logical(cat(1,Floe.alive)));
%plot_Floes_poly(0,0, Floe, ocean, c2_boundary);

%Calc high and low-res Eulerian fields
%[x,y, cFine0, cCoarse0,  U_Fine0,V_Fine0, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact );
%[c,vel,accel] = calc_eulerian_data(Floe,20,20,c2_boundary);
coarseMean=zeros(9,Ny,Nx,nSnapshots);
Vd = zeros(Ny,Nx,2);


%% Initialize time and other stuff to zero
if isempty(dir('dataDemAtmQG')); disp('Creating folder: dataDemAtmQG'); mkdir('dataDemAtmQG'); end
if isempty(dir('dataDemAtmQGc')); disp('Creating folder: dataDemAtmQGc'); mkdir('dataDemAtmQGc'); end

if ~exist('Time','var')
    Time=0;
    i_step=0;
    im_num=1;
    fig=0;
    %    EulCoarse=zeros(3, length(cCoarse0(:)),nSnapshots); %allocating memory
end

%psilowmin = -3; psilowmax = 3;

%% Solving for floe trajectories
tic;
gridArea=area(c2_boundary_poly)/Nx/Ny;
Vdnew=zeros(Ny, Nx);
figg = figure;

[psi12,pv12,qp] = AtmQG(0.7,pvold,i_step,1000,N,dt*nDTOut);
psi1 = real(ifft2(psi12));
psi_winds = psi1(:,:,2); % the second layer winds enteracting with ices
winds = UpdateWinds(psi_winds,Lx,Ly,dXa);
pvold = pv12;
        
ind = 1;
while im_num<nSnapshots
    
    display(i_step);
    
    if mod(i_step,nDTOut)==0  %plot the state after a number of timesteps
        
        % calculate winds for one step time-marching
        disp('Running atmosphere QG equations...');
        [psi12,pv12,qp] = AtmQG(0.7,pvold,i_step,10,N,dt*nDTOut);
        psi1 = real(ifft2(psi12));
        psi_winds = 0.5* ( psi1(:,:,1) + psi1(:,:,2) ); % barotropic mode
        %psi_winds = psi1(:,:,2); % the second layer winds enteracting with ices
        winds = UpdateWinds(psi_winds,Lx,Ly,dXa);
        pvold = pv12;
        
        ow = ocean;
        ow.Uocn = zeros(N+1); ow.Vocn = zeros(N+1);  % replace ocean by ow;
        ow.Uocn(1:N,1:N) = winds.U; ow.Vocn(1:N,1:N) = winds.V;
        
        figure(figg);
        figg=plot_Floes_poly_doublePeriodicBC(figg,Time,Floe, ow, c2_boundary_poly, PERIODIC); % plots model state
        saveas(figg,['./dataDemAtmQG/' num2str(im_num,'%03.f') '.jpg'],'jpg');
        
        %figure(figg);
        %figg=plotAtmWinds(figg,ow,Time);
        %saveas(figg,['./dataDemAtmQG/' num2str(im_num,'w%03.f') '.jpg'],'jpg');
        
        im_num = im_num+1;  %image number for saving data and coarse vars;
    end
    
    %Calculate forces and torques and intergrate forward
    [Floe,dissolvedNEW, SackedOB] = floe_interactions_all_doublePeriodicBCs_bpm(Floe, ocean, winds,heat_flux,c2_boundary_poly, dt,dissolvedNEW,SackedOB,Nx,Ny, RIDGING, PERIODIC,SUBFLOES);
    
    if FRACTURES
        overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
        keep=rand(length(Floe),1)>overlapArea;
        fracturedFloes=fracture_floe(Floe(~keep),3);
        %if length(fracturedFloes)<length(Floe(~keep)), disp('fractures killed floes'); end
        if ~isempty(fracturedFloes), fracturedFloes=rmfield(fracturedFloes,'potentialInteractions');
            Floe=[Floe(keep) fracturedFloes];
            %figure; plot([fracturedFloes.poly]); drawnow;
        end
    end
    
    %diluted=length(keep)-sum(keep);
    %if diluted>0, disp(['diluted floes: ' num2str(diluted)]); end
    
    [eularian_data] = calc_eulerian_data2(Floe,Nx,Ny,c2_boundary,PERIODIC);
    
    coarseMean(1,:,:,im_num)=squeeze(coarseMean(1,:,:,im_num))+eularian_data.c/nDTOut;
    coarseMean(2,:,:,im_num)=squeeze(coarseMean(2,:,:,im_num))+eularian_data.u/nDTOut;
    coarseMean(3,:,:,im_num)=squeeze(coarseMean(3,:,:,im_num))+eularian_data.v/nDTOut;
    coarseMean(4,:,:,im_num)=squeeze(coarseMean(4,:,:,im_num))+eularian_data.du/nDTOut;
    coarseMean(5,:,:,im_num)=squeeze(coarseMean(5,:,:,im_num))+eularian_data.dv/nDTOut;
    coarseMean(6,:,:,im_num)=squeeze(coarseMean(2,:,:,im_num))+eularian_data.mom_x/nDTOut;
    coarseMean(7,:,:,im_num)=squeeze(coarseMean(3,:,:,im_num))+eularian_data.mom_y/nDTOut;
    coarseMean(8,:,:,im_num)=squeeze(coarseMean(4,:,:,im_num))+eularian_data.force_x/nDTOut;
    coarseMean(9,:,:,im_num)=squeeze(coarseMean(5,:,:,im_num))+eularian_data.force_y/nDTOut;
    c = eularian_data.c;
    save(['./dataDemAtmQGc/c' num2str(ind,'%04.f') '.mat'],'c');  ind = ind+1;
    
    Area=cat(1,Floe.area);
    dissolvedNEW = calc_vol_dissolved(Floe(Area<3e5),Nx,Ny,c2_boundary_poly)+dissolvedNEW;
    if dhdt > 0
        dissolvedNEW = dissolvedNEW + (1-eularian_data.c)*gridArea*dhdt/(24*3600)*dt;
    end
    %Vd(:,:,im_num) = Vd(:,:,im_num)+Dissolved_Ice(Vd,coarseMean,im_num,dissolvedNEW,c2_boundary,dt)/nDTOut;
    Vdnew = Dissolved_Ice(Vd,coarseMean,im_num,dissolvedNEW,c2_boundary,dt);
    Vd(:,:,2) = Vd(:,:,1);
    Vd(:,:,1) = Vdnew;
    if PACKING
        [floe2,Vd] = pack_ice(Floe,c2_boundary,dhdt,Vd,SUBFLOES,PACKING);
    end
    %Vd(:,:,1)= Vd(:,:,1)+dissolvedNEW;
    
    Floe=Floe(Area> 1e6);
    if sum(Area<1e6)>0, display(['num of small floes killed:' num2str(sum(Area<1e6))]); end
    Time=Time+dt; i_step=i_step+1; %update time index
    
    for ii = 1:length(Floe)
        if Floe(ii).poly.NumRegions > 1
            xx(1) = 1;
            xx(1) = [1 2];
        end
    end
end
