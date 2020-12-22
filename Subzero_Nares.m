close all; clear all;

%% Set Flags

RIDGING=false; 

FRACTURES=true;

PERIODIC=true;

PACKING = false;

WELDING = false;

CORNERS = true;

ifPlot = false; %Plot floe figures or not?

%% Initialize model vars

dt=30; %Time step in sec
h0 = 0.1; %thickness of ice that gets packed in

%Define ocean currents
[ocean, HFo, nDTpack]=initialize_ocean_Nares(dt,h0);
nDTpack = 50;

%Define 10m winds
winds=[0 -10];

%Define boundaries
c2_boundary=initialize_boundaries();
Ly = max(c2_boundary(2,:));
c2_boundary_poly = polyshape(c2_boundary');
min_floe_size = 1e7;

%Initialize Floe state
%Floe=initialize_Floe('FloeShapes.mat');
height.mean = 2;
height.delta = 0;
target_concentration = [1;1;0];
[Floe, Nb] = initial_concentration_Nares(c2_boundary,target_concentration,height,50,min_floe_size);
% if isfield(Floe,'poly')
%     Floe=rmfield(Floe,{'poly'});
% end
%load Floe0000076.mat; Nb = 2;
%save('Floe0.mat','Floe')
%load('PackedFloesFullDomain.mat','Floe');
%Floe= create_packed_domain();

%%

dhdt = 1;

nDTOut=100; %Output frequency (in number of time steps)

nSnapshots=1000; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nPar = 4; %Number of workerks for parfor

target_concentration=1;

% specify coarse grid size
LxO= 2*max(ocean.Xo);LyO= 2*max(ocean.Yo);
Nx=10; Ny=fix(Nx*LyO/LxO);

%initialize dissolved ice at zero
dissolvedNEW=zeros(Ny,Nx);
Vd = zeros(Ny,Nx,2);
Vdnew=zeros(Ny, Nx);

%% Calc interactions and plot initial state
Floe=Floe(logical(cat(1,Floe.alive)));
[Floe,dissolvedNEW] = floe_interactions_all(Floe, ocean, winds,c2_boundary, dt,HFo,min_floe_size,Nx,Ny,Nb, dissolvedNEW,PERIODIC, RIDGING); % find interaction points
A=cat(1,Floe.area);
Amax = max(A);

%% Initialize time and other stuff to zero
if isempty(dir('figs')); disp('Creating folder: figs'); mkdir('figs'); end
if isempty(dir('FloesNares')); disp('Creating folder: Floes'); mkdir('FloesNares'); end

if ~exist('Time','var')
    Time=0;
    i_step=0;
    im_num=1;
    fig=0;
end


%% Solving for floe trajectories
Floe3 = Floe;
Floe2 = Floe;
tic;
while im_num<nSnapshots

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

        if ifPlot
            [fig] =plot_basic(fig, Time,Floe,ocean,c2_boundary_poly,Nb);
            saveas(fig,['./figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
        end
        
        save(['./FloesNares/Floe' num2str(im_num,'%07.f') '.mat'],'Floe');
        
        im_num=im_num+1;  %image number for saving data and coarse vars;
    end

    if PACKING && h0 > 0
        if mod(i_step,nDTpack)==0corn
            height.mean = h0;
            height.delta = 0;
            [Floe,Vd] = pack_ice(Floe,c2_boundary,dhdt,Vd,target_concentration,ocean,height, min_floe_size, PERIODIC);
        end
    end
    
%     Floe3 = Floe2; Floe2 = Floe;
    %Calculate forces and torques and intergrate forward
    [Floe,dissolvedNEW] = floe_interactions_all(Floe, ocean, winds, c2_boundary, dt, HFo,min_floe_size, Nx,Ny,Nb, dissolvedNEW,PERIODIC, RIDGING);
    
    if FRACTURES
        overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
        keep=rand(length(Floe),1)>overlapArea;
        keep(1:Nb) = ones(Nb,1);
        fracturedFloes=fracture_floe(Floe(~keep),3);
        if ~isempty(fracturedFloes)
            Floe=[Floe(keep) fracturedFloes];
        end
        
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
        keep=stresses<rand(length(Floe),1);
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
%%
