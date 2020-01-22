close all; clear all;


%% Initialize model vars

%Define ocean currents
ocean=initialize_ocean();

%Define 10m winds
winds=[0 5];

%Initialize Floe state
%Floe=initialize_Floe('FloeShapes.mat');
%load('Floe.mat','Floe');
load('PackedFloesFullDomain.mat','Floe');

Floe= create_polygons(Floe);
Floe=RemoveFullyOverlappingFloes(Floe);
Floe=rmfield(Floe,{'interactions','potentialInteractions'});
Floe=rmfield(Floe,{'c0','c_alpha','X','Y'});


%Floe= create_packed_domain();
%Define boundaries
c2_boundary=initialize_boundaries();
in=inpolygon(cat(1,Floe.Xi),cat(1,Floe.Yi),c2_boundary(1,:),c2_boundary(2,:));
Floe=Floe(logical(in));
%%

dt=40; %Time step in sec

nDTOut=50; %Output frequency (in number of time steps)

nSnapshots=1000; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nPar = 4; %Number of workerks for parfor

ifPlot = true; %Plot floe figures or not?


%% Calc interactions and plot initial state
Floe = floe_interactions_all(Floe, ocean, winds,c2_boundary, dt); % find interaction points
Floe=Floe(logical(cat(1,Floe.alive)));
plot_Floes_poly(0,0, Floe, ocean, c2_boundary);

%% Define Eulerian grid and coarsening factor
ddx=250; % resolution of the original floe images in meters
[Xgg, Ygg]=meshgrid(-70e3:ddx:70e3,-70e3:ddx:70e3); % high-res eulerian grid
c_fact=40; % coarsening factor

%Calc high and low-res Eulerian fields
%[x,y, cFine0, cCoarse0,  U_Fine0,V_Fine0, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact );


%% Initialize time and other stuff to zero
if isempty(dir('figs')); disp('Creating folder: figs'); mkdir('figs'); end

if ~exist('Time','var')
    Time=0;
    i_step=0;
    im_num=1;
    fig=0;
%    EulCoarse=zeros(3, length(cCoarse0(:)),nSnapshots); %allocating memory
end


%% Solving for floe trajectories
tic;
while im_num<nSnapshots
     
    %c2_boundary=c2_boundary*(1+0.0005); % shrink by % every 10 steps
    display(i_step);
    if mod(i_step,10)==0
        disp(newline);
        toc
        disp([num2str(i_step) ' timesteps comleted']); 
        numCollisions = calc_collisionNum(Floe);
        sacked = sum(~cat(1, Floe.alive));
        if sacked>0, disp(['sacked floes: ' num2str(sacked)]); end
        disp(['number of collisions: ' num2str(numCollisions)  newline]);
        tic
    end

    if mod(i_step,nDTOut)==0  %plot the state after a number of timesteps
        
        Floe=RemoveFullyOverlappingFloes(Floe);  % remove highly overlapping elements

        
        if ifPlot
            fig=plot_Floes_poly(fig,Time,Floe, ocean, c2_boundary); % plots model state
            saveas(fig,['./figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
        end
        
        %calculating and saving corase grid variables
        %[x,y, cFine0, cCoarse0,  U_Fine0,V_Fine0, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact );
%         [~,~, ~, cCoarse0,  ~,~, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact );
%         EulCoarse(1,:,im_num)= cCoarse0(:);
%         EulCoarse(2,:,im_num)= U_Coarse0(:);
%         EulCoarse(3,:,im_num)= V_Coarse0(:);
        
        im_num=im_num+1;  %image number for saving data and coarse vars;
    end
    
    
    %Calculate forces and torques and intergrate forward
    Floe = floe_interactions_all(Floe, ocean, winds, c2_boundary, dt);
    
    Time=Time+dt; i_step=i_step+1; %update time index


end
%%
