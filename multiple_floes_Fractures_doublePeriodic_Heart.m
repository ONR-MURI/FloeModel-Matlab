close all; clear all;

addpath ~/Downloads/dengwirda-inpoly-ebf47d6/ 

%% Initialize model vars

%Define ocean currents
[ocean, c2_boundary]=initialize_ocean_Gyre(1e4, 2e5, 1e5,4e3);
c2_boundary_poly=polyshape(c2_boundary(1,:),c2_boundary(2,:));

c=1; % could be a vector
Floe = initialize_concentration(c,c2_boundary,100);


t = linspace(-pi,pi, 30);
X =      t .* sin( pi * .872*sin(t)./t);
Y = -abs(t) .* cos(pi * sin(t)./t);

c2_boundary_heart=polyshape(5e3*X,5e3*(Y+1));

for i=1:length(Floe)

poly=subtract(Floe(i).poly,c2_boundary_heart);
poly=regions(poly);
if length(poly)==1
    Floe(i)=rmfield(initialize_floe_values(poly),'potentialInteractions');
else
    Floe(i)=[];
    for k=1:length(poly)
    floe=rmfield(initialize_floe_values(poly(k)),'potentialInteractions');
    Floe=[Floe floe];
    end
end

end

floe0=rmfield(initialize_floe_values(c2_boundary_heart),'potentialInteractions');
Floe=[floe0 Floe];
%Floe=Floe(~isempty(Floe));
%Define 10m winds
winds=[0 0];

%Initialize Floe state
%load('Floe_clean.mat','Floe');


plot_Floes_poly(0,0, Floe, ocean, c2_boundary);
%%

dt=20; %Time step in sec

nDTOut=10; %Output frequency (in number of time steps)

nSnapshots=10000; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nPar = 4; %Number of workerks for parfor

ifPlot = true; %Plot floe figures or not?


% Calc interactions and plot initial state
Floe = floe_interactions_all_doublePeriodicBCs_Heart(Floe, ocean, winds,c2_boundary_poly, dt); % find interaction points
Floe=Floe(logical(cat(1,Floe.alive)));
%plot_Floes_poly(0,0, Floe, ocean, c2_boundary);

%% Define Eulerian grid and coarsening factor
%ddx=250; % resolution of the original floe images in meters
%[Xgg, Ygg]=meshgrid(-70e3:ddx:70e3,-70e3:ddx:70e3); % high-res eulerian grid
%c_fact=40; % coarsening factor

%Calc high and low-res Eulerian fields
%[x,y, cFine0, cCoarse0,  U_Fine0,V_Fine0, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact );
[c,vel,accel] = calc_eulerian_data(Floe,20,20,c2_boundary);

Nx=10; Ny=10;

coarseMean=zeros(5,Ny,Nx,nSnapshots);
coarseSnap=zeros(5,Ny,Nx,nSnapshots);
Vd = zeros(Ny,Nx,2);

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
gridArea=area(c2_boundary_poly)/Nx/Ny;
Vdnew=zeros(Ny,Nx);
fig2=figure;
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
                
        if ifPlot
            fig=plot_Floes_poly_doublePeriodicBC(fig,Time,Floe, ocean, c2_boundary_poly); % plots model state
            saveas(fig,['./figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
            if im_num>1
            if (~isvalid(fig2)), fig2=figure; end
            figure(fig2);
            imagesc(Vdnew/gridArea/1e3); axis xy
            u=squeeze(coarseMean(2,:,:,im_num));
            v=squeeze(coarseMean(3,:,:,im_num));
            colormap('gray');colorbar;
            hold on; quiver(u,v,'r')
            drawnow
            figure(fig);
            end

        end
        
        %calculating and saving corase grid variables
        
        [c,vel,accel] = calc_eulerian_data(Floe,Nx,Ny,c2_boundary);
        coarseSnap(1,:,:,im_num)=c;
        coarseSnap(2,:,:,im_num)=vel.u;
        coarseSnap(3,:,:,im_num)=vel.v;
        coarseSnap(4,:,:,im_num)=accel.du;
        coarseSnap(5,:,:,im_num)=accel.dv;
        
        save('coarseData.mat','coarseSnap','coarseMean');
        
        %calculating and saving corase grid variables
        %[x,y, cFine0, cCoarse0,  U_Fine0,V_Fine0, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact );
%         [~,~, ~, cCoarse0,  ~,~, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact );
%         EulCoarse(1,:,im_num)= cCoarse0(:);
%         EulCoarse(2,:,im_num)= U_Coarse0(:);
%         EulCoarse(3,:,im_num)= V_Coarse0(:);
        
        im_num=im_num+1;  %image number for saving data and coarse vars;
    end
    
    %Calculate forces and torques and intergrate forward
    Floe = floe_interactions_all_doublePeriodicBCs_Heart(Floe, ocean, winds, c2_boundary_poly, dt);
    %Floe = floe_interactions_all_periodicBCs_bpm(Floe, ocean, winds, c2_boundary_poly, dt,dissolvedNEW,Nx,Ny);

    overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area); 
    overlapArea(1)=0; %keeping floe#1 intact
    keep=rand(length(Floe),1)>overlapArea;
    fracturedFloes=fracture_floe(Floe(~keep),3);        
    %if length(fracturedFloes)<length(Floe(~keep)), disp('fractures killed floes'); end
    if ~isempty(fracturedFloes), fracturedFloes=rmfield(fracturedFloes,'potentialInteractions'); 
    Floe=[Floe(keep) fracturedFloes];
    %figure; plot([fracturedFloes.poly]); drawnow;
    end
    %diluted=length(keep)-sum(keep);
    %if diluted>0, disp(['diluted floes: ' num2str(diluted)]); end
    
    [c,vel,accel] = calc_eulerian_data(Floe,Nx,Ny,c2_boundary);
    
    coarseMean(1,:,:,im_num)=squeeze(coarseMean(1,:,:,im_num))+c/nDTOut;
    coarseMean(2,:,:,im_num)=squeeze(coarseMean(2,:,:,im_num))+vel.u/nDTOut;
    coarseMean(3,:,:,im_num)=squeeze(coarseMean(3,:,:,im_num))+vel.v/nDTOut;
    coarseMean(4,:,:,im_num)=squeeze(coarseMean(4,:,:,im_num))+accel.du/nDTOut;
    coarseMean(5,:,:,im_num)=squeeze(coarseMean(5,:,:,im_num))+accel.dv/nDTOut;
    
    Area=cat(1,Floe.area);
    dissolvedNEW = calc_vol_dissolved(Floe(Area<3e5),Nx,Ny,c2_boundary_poly);
    %Vd(:,:,im_num) = Vd(:,:,im_num)+Dissolved_Ice(Vd,coarseMean,im_num,dissolvedNEW,c2_boundary,dt)/nDTOut;
    Vdnew = Dissolved_Ice(Vd,coarseMean,im_num,dissolvedNEW,c2_boundary,dt);
    Vd(:,:,2) = Vd(:,:,1);
   % Vd(:,:,1) = Vdnew;
    Vd(:,:,1)= Vd(:,:,1)+dissolvedNEW;
    
    Floe=Floe(Area> 1e6);
    if sum(Area<1e6)>0, display(['num of small floes killed:' num2str(sum(Area<1e6))]); end
    Time=Time+dt; i_step=i_step+1; %update time index


end
%%
