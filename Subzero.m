close all; clear all;

%% Set Flags

RIDGING=false; 

FRACTURES=false;

PERIODIC=false;

PACKING = false;

WELDING = false;

ifPlot = false; %Plot floe figures or not?

%% Initialize model vars

%Define ocean currents
ocean=initialize_ocean();

%Define 10m winds
winds=[10 10];

%Define boundaries
c2_boundary=initialize_boundaries();
Ly = max(c2_boundary(2,:));
c2_boundary_poly = polyshape(c2_boundary');
Nb = 0;
min_floe_size = 1e6;

%Initialize Floe state
%Floe=initialize_Floe('FloeShapes.mat');
height.mean = 2;
height.delta = 0.25;
target_concentration = 0.75;
Floe = initial_concentration(c2_boundary,target_concentration,height,100,min_floe_size);
if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
end
%load('PackedFloesFullDomain.mat','Floe');
%Floe= create_packed_domain();

%%

dt=10; %Time step in sec

h0 = 0.1;

dhdt = 1;

nDTOut=50; %Output frequency (in number of time steps)

nSnapshots=1000; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nPar = 4; %Number of workerks for parfor

target_concentration=1;

Vd = zeros(5,5,2);

%% Calc interactions and plot initial state
Floe=Floe(logical(cat(1,Floe.alive)));
for ii = 1:length(Floe)
    Floe(ii).h = 2;
end
Floe = floe_interactions_all(Floe, ocean, winds,c2_boundary, dt,PERIODIC, RIDGING); % find interaction points
Nb = 0;
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
end


%% Solving for floe trajectories
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
        
        if ifPlot
            [fig] =plot_basic(fig, Time,Floe,ocean,c2_boundary_poly,Nb);
            saveas(fig,['./figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
        end
        
        simp = 0;
        for ii = 1:length(Floe)
            if length(Floe(ii).c0) > 30
                Floe(ii) = FloeSimplify(Floe(ii));
                simp = simp+1;
            end
        end
        
        if PACKING && h0 > 0
            height.mean = h0;
            height.delta = 0;
            [Floe,Vd] = pack_ice(Floe,c2_boundary,dhdt,Vd,target_concentration,ocean,height, min_floe_size, PERIODIC);
        end
        
        if WELDING
            weldrate = 0.1;%Set rate at which floes will meld
            A=cat(1,Floe.area);
            if max(A) > Amax
                Amax = max(A);
            end
            FloeOld = Floe;
            Floe = Weld_Floes(Floe,weldrate,Amax);
        end
        
        save(['./Floes/Floe' num2str(im_num,'%07.f') '.mat'],'Floe');
        
        im_num=im_num+1;  %image number for saving data and coarse vars;
    end
    
    
    %Calculate forces and torques and intergrate forward
    Floe = floe_interactions_all(Floe, ocean, winds, c2_boundary, dt, PERIODIC, RIDGING);
    
    if FRACTURES
        overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
        keep=rand(length(Floe),1)>2*overlapArea;
        fracturedFloes=fracture_floe(Floe(~keep),3);
        if ~isempty(fracturedFloes)
            Floe=[Floe(keep) fracturedFloes];
        end
        
    end
    
%     Floe = corners(Floe);
    
    Area=cat(1,Floe.area);
    if sum(Area<min_floe_size)>0, display(['num of small floes killed:' num2str(sum(Area<min_floe_size))]); end
    Floe=Floe(Area> min_floe_size);
    live = cat(1,Floe.alive);
    Floe(live == 0) = [];
    
    Time=Time+dt; i_step=i_step+1; %update time index

end
%%
