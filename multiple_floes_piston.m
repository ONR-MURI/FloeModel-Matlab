close all; clear all;


%% Initialize model vars

%Define ocean currents
ocean=initialize_ocean_piston();

%Define 10m winds
winds=[0 0];

%Initialize Floe state
%Floe=initialize_Floe('FloeShapes.mat');
load('FloeVoronoi.mat','Floe');
%load('PackedFloesFullDomain.mat','Floe');
%Floe= create_packed_domain();
%Define boundaries
c2_boundary=initialize_boundaries();
Ly = max(c2_boundary(2,:));
c2_boundary_poly = polyshape(c2_boundary');
Nb = 0;

%%

dt=10; %Time step in sec

nDTOut=50; %Output frequency (in number of time steps)

nSnapshots=1000; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nPar = 4; %Number of workerks for parfor

ifPlot = false; %Plot floe figures or not?

min_floe_size = 1e6;

%% Calc interactions and plot initial state
Floe=Floe(logical(cat(1,Floe.alive)));
for ii = 1:length(Floe)
    Floe(ii).h = 2;
end
Floe = floe_interactions_all(Floe, ocean, winds,c2_boundary, dt); % find interaction points
%plot_Floes(0,0, Floe, ocean, c2_boundary);


%% Define Eulerian grid and coarsening factor
ddx=250; % resolution of the original floe images in meters
[Xgg, Ygg]=meshgrid(-70e3:ddx:70e3,-70e3:ddx:70e3); % high-res eulerian grid
c_fact=40; % coarsening factor

%Calc high and low-res Eulerian fields
[x,y, cFine0, cCoarse0,  U_Fine0,V_Fine0, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact,c2_boundary );




%% Initialize time and other stuff to zero
if isempty(dir('figs')); disp('Creating folder: figs'); mkdir('figs'); end
if isempty(dir('Floes')); disp('Creating folder: Floes'); mkdir('Floes'); end

if ~exist('Time','var')
    Time=0;
    i_step=0;
    im_num=1;
    fig=0;
    EulCoarse=zeros(3, length(cCoarse0(:)),nSnapshots); %allocating memory
end


%% Solving for floe trajectories
tic;
while im_num<nSnapshots
     
     
    Ly = Ly-5; % shrink the boundaries to crush from top and bottom
    c2_boundary(2,:) = [-Ly Ly Ly -Ly -Ly];% 
    
%     c2_boundary(1,:) = c2_boundary(1,:) + [-10 10 10 -10 -10]; %Shearing boundary

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
            %fig=plot_Floes(fig,Time,Floe, ocean, c2_boundary); % plots model state
            saveas(fig,['./figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
        end
        
        simp = 0;
        for ii = 1:length(Floe)
            if length(Floe(ii).c0) > 50
                Floe(ii) = FloeSimplify(Floe(ii));
                simp = simp+1;
            end
        end
        save(['./Floes/Floe' num2str(im_num,'%07.f') '.mat'],'Floe');

        %calculating and saving corase grid variables
        %[x,y, cFine0, cCoarse0,  U_Fine0,V_Fine0, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact,c2_boundary );
%         [~,~, ~, cCoarse0,  ~,~, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact ,c2_boundary);
%         EulCoarse(1,:,im_num)= cCoarse0(:);
%         EulCoarse(2,:,im_num)= U_Coarse0(:);
%         EulCoarse(3,:,im_num)= V_Coarse0(:);
        
        im_num=im_num+1;  %image number for saving data and coarse vars;
    end
    
    
    %Calculate forces and torques and intergrate forward
    Floe = floe_interactions_all(Floe, ocean, winds, c2_boundary, dt);
    
    FRACTURES = true;
    if FRACTURES
        overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
        keep=rand(length(Floe),1)>2*overlapArea;
        fracturedFloes=fracture_floe(Floe(~keep),3);
        if ~isempty(fracturedFloes)
            Floe=[Floe(keep) fracturedFloes];
        end
        
    end
    
    Area=cat(1,Floe.area);
    if sum(Area<min_floe_size)>0, display(['num of small floes killed:' num2str(sum(Area<min_floe_size))]); end
    Floe=Floe(Area> min_floe_size);
    
    Time=Time+dt; i_step=i_step+1; %update time index

end
%%
