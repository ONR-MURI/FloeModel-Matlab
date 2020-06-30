close all; clear all;

addpath ~/Downloads/dengwirda-inpoly-ebf47d6/ 

%% Initialize model vars
RIDGING=false; 

FRACTURES=false;

PERIODIC=true;

SUBFLOES = false;

PACKING = false;

WELDING = false;

%Define ocean currents
transport=0*1e4;
Lx=2e5; Ly=1e5; dXo=4e3;
[ocean, c2_boundary]=initialize_ocean_Gyre(transport, Lx, Ly,dXo);
c2_boundary_poly=polyshape(c2_boundary(1,:),c2_boundary(2,:));

%Define 10m winds
winds=[5 0];

%Initialize Floe state
height.mean = 2;
height.delta = 1; %max difference between a flow thickness and the mean floe value

target_concentration=0.7; % could be a vector
Floe = initialize_concentration(target_concentration, c2_boundary,ocean,SUBFLOES,height, 75);
% load FloeBase
%plot_Floes_poly(0,0, Floe, ocean, c2_boundary);
%%

dt=20; %Time step in sec

nDTOut=10; %Output frequency (in number of time steps)

nSnapshots=10000; %Total number of model snapshots to save

nDT=nDTOut*nSnapshots; %Total number of time steps

nPar = 4; %Number of workerks for parfor

ifPlot = true; %Plot floe figures or not?

SackedOB = 0; %initialize number of floes sacked for being out of bounds at zero

dhdt = 0; %Sets probability for ice replenishing open space

heat_flux = -2.1*18/(334*1000*1000)*(100*3600*24); %Rate at which ice thickness increases thermodynamically in cm/day (once divided by h)
heat_flux = 0;

% specify coarse grid size
Lx= 2*max(ocean.Xo);Ly= 2*max(ocean.Yo);
Nx=20; Ny=fix(Nx*Ly/Lx);

%Track floe states
%NumFloes = zeros(1,nSnapshots);
Areas = zeros(2,10,nSnapshots);
Thicknesses = zeros(2,10,nSnapshots);
NumFloes = length(Floe);
Atot = repmat(cat(1,Floe.area),nDTOut);
htot = repmat(cat(1,Floe.h),nDTOut);

%initialize dissolved ice at zero
dissolvedNEW=zeros(Ny,Nx);
Vd = zeros(Ny,Nx,2);

% Calc interactions and plot initial state
[Floe,dissolvedNEW, SackedOB] = floe_interactions_all_doublePeriodicBCs_bpm(Floe, ocean, winds,heat_flux,c2_boundary_poly, dt,dissolvedNEW,SackedOB,Nx,Ny, RIDGING, PERIODIC,SUBFLOES); % find interaction points
Floe=Floe(logical(cat(1,Floe.alive)));
%plot_Floes_poly(0,0, Floe, ocean, c2_boundary);

%% Define Eulerian grid and coarsening factor
%ddx=250; % resolution of the original floe images in meters
%[Xgg, Ygg]=meshgrid(-70e3:ddx:70e3,-70e3:ddx:70e3); % high-res eulerian grid
%c_fact=40; % coarsening factor

%Calc high and low-res Eulerian fields
%[x,y, cFine0, cCoarse0,  U_Fine0,V_Fine0, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact );
[eularian_data] = calc_eulerian_data2(Floe,Nx,Ny,c2_boundary,PERIODIC);

coarseMean=zeros(9,Ny,Nx,nSnapshots);
coarseSnap=zeros(9,Ny,Nx,nSnapshots);
A=cat(1,Floe.area);
Amax = max(A);
SimpMin = @(A) 15*log10(A);%

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
Vdnew=zeros(Ny, Nx);
fig2=figure;
fig3=figure;
while im_num<nSnapshots
     
    %c2_boundary=c2_boundary*(1+0.0005); % shrink by % every 10 steps
    %display(i_step);
    if mod(i_step,10)==0
        disp(newline);
        toc
        disp([num2str(i_step) ' timesteps comleted']); 
        numCollisions = calc_collisionNum(Floe);
        sacked = sum(~cat(1, Floe.alive));
        if sacked>0, disp(['sacked floes: ' num2str(sacked)]); end
        if SackedOB>0, disp(['total sacked floes for being out of bounds: ' num2str(SackedOB)]); end
        disp(['number of collisions: ' num2str(numCollisions)  newline]);
        tic
    end

    if mod(i_step,nDTOut)==0  %plot the state after a number of timesteps
                
        
        
        %calculating and saving corase grid variables
        
        [eularian_data] = calc_eulerian_data2(Floe,Nx,Ny,c2_boundary,PERIODIC);
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
        save('Floe.mat','Floe');
        
        %NumFloes(im_num) = length(Floe);
        %A = cat(1,Floe.area);
        [A1,A2] = hist(Atot,10);
        %h = cat(1,Floe.h);
        [h1,h2] = hist(htot,10);
        FloeStats(im_num).Num = fix(NumFloes);
        FloeStats(im_num).A = Atot;
        FloeStats(im_num).h = htot;
        %Areas(:,:,im_num) = [A1/nDTOut;A2];
        %Thicknesses(:,:,im_num) = [h1/nDTOut;h2];
        Atot = [];
        htot = [];
        NumFloes = 0;
        save('FloeStats.mat','FloeStats','Amax')
        
        if mod(i_step,nDTOut)==0 && PACKING
            height.mean = 0.2;
            height.delta = 0;
            [Floe,Vd] = pack_ice(Floe,c2_boundary,dhdt,Vd,target_concentration,ocean,SUBFLOES, PERIODIC);
        end

        floenew = [];
        for ii = 1:length(Floe)
            floe = Floe(ii);
            ddx = 100;
            while length(floe(1).poly.Vertices) > SimpMin(Floe(ii).area)
                floe = FloeSimplify(Floe(ii), ddx,SUBFLOES);
                ddx = ddx + 150;
            end
            for jj = 1:length(floe)
                if jj == 1
                    Floe(ii) = floe(jj);
                else
                    floenew = [floenew floe(jj)];
                end
            end           
        end
        Floe = [Floe floenew];
        A = cat(1,Floe.area);
        Floe(A<3500) = [];
        live = cat(1,Floe.alive);
        Floe(live==0)=[];
        
        if ifPlot
            fig=plot_Floes_poly_doublePeriodicBC(fig,Time,Floe, ocean, c2_boundary_poly, PERIODIC); % plots model state
            %saveas(fig,['./figs/' num2str(im_num,'%03.f') '.jpg'],'jpg');
            %figure(fig3);
            %fig3=plot_Floes_poly_doublePeriodicBC_thickness(fig3,Time,Floe, ocean, c2_boundary_poly, PERIODIC); 
            %saveas(fig3,['./figs/' num2str(im_num,'t%03.f') '.jpg'],'jpg');
            if im_num>1
            if (~isvalid(fig2)), fig2=figure; end
            figure(fig2);
            %imagesc(Vdnew/gridArea/1e3); axis xy
            imagesc(1:Nx,1:Ny,eularian_data.c); axis ij
            colormap('gray'); colorbar;
            hold on; quiver(eularian_data.u,eularian_data.v,'r')
            drawnow
            end

        end
        
        %calculating and saving corase grid variables
        %[x,y, cFine0, cCoarse0,  U_Fine0,V_Fine0, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact );
%         [~,~, ~, cCoarse0,  ~,~, U_Coarse0, V_Coarse0 ] = create_eulerian_data( Floe, Xgg, Ygg, c_fact );
%         EulCoarse(1,:,im_num)= cCoarse0(:);
%         EulCoarse(2,:,im_num)= U_Coarse0(:);
%         EulCoarse(3,:,im_num)= V_Coarse0(:);
        
        im_num=im_num+1;  %image number for saving data and coarse vars;
    end
    
    save('FloeOld.mat','Floe')
    %Calculate forces and torques and intergrate forward
    [Floe,dissolvedNEW, SackedOB] = floe_interactions_all_doublePeriodicBCs_bpm(Floe, ocean, winds,heat_flux, c2_boundary_poly, dt,dissolvedNEW,SackedOB,Nx,Ny, RIDGING, PERIODIC,SUBFLOES);
    
%     save('FloeNow2.mat','Floe')
    
    if mod(i_step-1,nDTOut)==0
        if WELDING
            meldrate = 0.2;%Set rate at which floes will meld
            A=cat(1,Floe.area);
            if max(A) > Amax
                Amax = max(A);
            end
            Floe = MeldFloes(Floe,meldrate,Amax,SUBFLOES,Nx,Ny,c2_boundary,PERIODIC);
        end
    end
    
    if FRACTURES
        overlapArea=cat(1,Floe.OverlapArea)./cat(1,Floe.area);
        keep=rand(length(Floe),1)>overlapArea;
        fracturedFloes=fracture_floe(Floe(~keep),5);
%         Floe=rmfield(Floe,{'potentialInteractions'});
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
    coarseMean(6,:,:,im_num)=squeeze(coarseMean(6,:,:,im_num))+eularian_data.mom_x/nDTOut;
    coarseMean(7,:,:,im_num)=squeeze(coarseMean(7,:,:,im_num))+eularian_data.mom_y/nDTOut;
    coarseMean(8,:,:,im_num)=squeeze(coarseMean(8,:,:,im_num))+eularian_data.force_x/nDTOut;
    coarseMean(9,:,:,im_num)=squeeze(coarseMean(9,:,:,im_num))+eularian_data.force_y/nDTOut;
    
    %NumFloes(im_num) =NumFloes(im_num)+ length(Floe)/nDTOut;
    NumFloes = NumFloes+length(Floe)/nDTOut;
    Asnap = cat(1,Floe.area);
    Atot = [Atot; Asnap];
    hsnap = cat(1,Floe.h);
    htot = [htot; hsnap];
    
    Area=cat(1,Floe.area);
    dissolvedNEW = calc_vol_dissolved(Floe(Area<3e5),Nx,Ny,c2_boundary_poly)+dissolvedNEW;
    if dhdt > 0
        dissolvedNEW = dissolvedNEW + (1-eularian_data.c)*gridArea*dhdt/(24*3600)*dt;
    end
    %Vd(:,:,im_num) = Vd(:,:,im_num)+Dissolved_Ice(Vd,coarseMean,im_num,dissolvedNEW,c2_boundary,dt)/nDTOut;
    Vdnew = Dissolved_Ice(Vd,coarseMean,im_num,dissolvedNEW,c2_boundary,dt);
    Vd(:,:,2) = Vd(:,:,1);
    Vd(:,:,1) = Vdnew;
    
     
    %Vd(:,:,1)= Vd(:,:,1)+dissolvedNEW;
    
    Floe=Floe(Area> 1e6);
    if sum(Area<1e6)>0, display(['num of small floes killed:' num2str(sum(Area<1e6))]); end
    Time=Time+dt; i_step=i_step+1; %update time index

end


