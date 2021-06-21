function [eularian_data] = calc_eulerian_stress(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC)
%% Function to take information of all floes and average them over a corase grained area
id = 'MATLAB:polyshape:boolOperationFailed';
warning('off',id)

%Identify only the live floes
live = cat(1,Floe.alive);
Floe(live==0)=[];
N0 = length(Floe);
for ii = 1:N0
    Stress(ii) = max(abs(eig(Floe(ii).Stress)));
    if isempty(Floe(ii).Fx)
        Floe(ii).Fx = 0;
        Floe(ii).Fy = 0;
    end
end
[~,TF] = rmoutliers(Stress);
for ii = 1:N0
    if TF(ii)
        Floe(TF(ii)).alive = 0;
    end
end
clear Stress

%Create ghost floes for periodic floe states
Lx= max(c2_boundary(1,:));
Ly= max(c2_boundary(2,:));%c2 must be symmetric around x=0 for channel boundary conditions.
N0=length(Floe);
if PERIODIC
    
    ghostFloeX=[];
    ghostFloeY=[];
    parent=[];
    translation = [];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        poly = polyshape(Floe(i).c_alpha'+[x(i) y(i)]);
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(poly.Vertices(:,1)))>Lx)
            
            ghostFloeX=[ghostFloeX  Floe(i)];
            ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
            
        end
        
        
    end
    
    Floe=[Floe ghostFloeX];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(poly.Vertices(:,2)))>Ly)
            
            ghostFloeY=[ghostFloeY  Floe(i)];
            ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
            
        end
        
    end
    
    Floe=[Floe ghostFloeY];
    
end

Floe(1:Nb) = [];

% Idenfity floes that are alive
live = cat(1,Floe.alive);

%Create coarse grid and coarse floe variables
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
y = fliplr(y);
dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
Area = dx*dy;
eularian_data.u = zeros(Ny,Nx);
eularian_data.v = zeros(Ny,Nx);
eularian_data.du = zeros(Ny,Nx);
eularian_data.dv = zeros(Ny,Nx);
eularian_data.force_x = zeros(Ny,Nx);
eularian_data.force_y = zeros(Ny,Nx);
eularian_data.stress = zeros(Ny,Nx);
eularian_data.stressxx = zeros(Ny,Nx);
eularian_data.stressyx = zeros(Ny,Nx);
eularian_data.stressxy = zeros(Ny,Nx);
eularian_data.stressyy = zeros(Ny,Nx);
eularian_data.strainux = zeros(Ny,Nx);
eularian_data.strainvx = zeros(Ny,Nx);
eularian_data.strainuy = zeros(Ny,Nx);
eularian_data.strainvy = zeros(Ny,Nx);


mass = cat(1,Floe.mass);
mass(isnan(mass)==1)=0;
A = cat(1,Floe.area);
A(isnan(A)==1)=0;
U = cat(1,Floe.Ui);
U(isnan(U)==1)=0;
V = cat(1,Floe.Vi);
V(isnan(V)==1)=0;
dU = cat(1,Floe.dUi_p);%(U-cat(1,Floe.dXi_p))/dt;
dU(isnan(dU)==1)=0;
dV = cat(1,Floe.dVi_p);%(V-cat(1,Floe.dYi_p))/dt;
dV(isnan(dV)==1)=0;
ForceX = cat(1,Floe.Fx);
ForceX(isnan(ForceX)==1)=0;
ForceY = cat(1,Floe.Fy);
ForceY(isnan(ForceY)==1)=0;
Stress = zeros(2,2,length(Floe));
Strain = zeros(2,2,length(Floe));
for ii = 1:length(Floe)
    Stress(:,:,ii) = Floe(ii).Stress;
    Strain(:,:,ii) = Floe(ii).strain;
end

%Find floes and create bins 
Xi=cat(1,Floe.Xi);
Yi=cat(1,Floe.Yi);
Binx = fix((Xi-min(x))/(max(x)-min(x))*Nx+1);
Biny = fix((Yi-min(y))/(max(y)-min(y))*Ny+1);

%% Loop to find coarse averages

for ii = 1:Nx
    for jj = 1:Ny
        Mtot = sum(cat(1,Floe(live == 1 & Binx == ii & Biny == jj).mass));
        if Mtot>0
            eularian_data.u(jj,ii) = sum(U(live == 1 & Binx == ii & Biny == jj).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;
            eularian_data.v(jj,ii) = sum(V(live == 1 & Binx == ii & Biny == jj).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;
            eularian_data.du(jj,ii) = sum(dU(live == 1 & Binx == ii & Biny == jj).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;
            eularian_data.dv(jj,ii) = sum(dV(live == 1 & Binx == ii & Biny == jj).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;
            eularian_data.force_x(jj,ii) = sum(ForceX(live == 1 & Binx == ii & Biny == jj))./Mtot;%sum(mass(logical(potentialInteractions(jj,ii,:)))'.*dU(logical(potentialInteractions(jj,ii,:)))'.*Aover)./(dt*Area);
            eularian_data.force_y(jj,ii) = sum(ForceY(live == 1 & Binx == ii & Biny == jj).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;%sum(mass(logical(potentialInteractions(jj,ii,:)))'.*dV(logical(potentialInteractions(jj,ii,:)))'.*Aover)./(dt*Area);
%             eularian_data.stress(jj,ii) = 0.5*trace(sum(Stress(:,:,logical(potentialInteractions(jj,ii,:))).*zg,3)./sum(Aover))/2;
            eularian_data.stressxx(jj,ii) = sum(squeeze(Stress(1,1,(live == 1 & Binx == ii & Biny == jj))).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;
            eularian_data.stressyx(jj,ii) = sum(squeeze(Stress(1,2,(live == 1 & Binx == ii & Biny == jj))).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;
            eularian_data.stressxy(jj,ii) = sum(squeeze(Stress(2,1,(live == 1 & Binx == ii & Biny == jj))).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;
            eularian_data.stressyy(jj,ii) = sum(squeeze(Stress(2,2,(live == 1 & Binx == ii & Biny == jj))).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;
            eularian_data.strainux(jj,ii) = sum(squeeze(Strain(1,1,(live == 1 & Binx == ii & Biny == jj))).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;
            eularian_data.strainvx(jj,ii) = sum(squeeze(Strain(1,2,(live == 1 & Binx == ii & Biny == jj))).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;
            eularian_data.strainuy(jj,ii) = sum(squeeze(Strain(2,1,(live == 1 & Binx == ii & Biny == jj))).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;
            eularian_data.strainvy(jj,ii) = sum(squeeze(Strain(2,2,(live == 1 & Binx == ii & Biny == jj))).*mass(live == 1 & Binx == ii & Biny == jj))./Mtot;
            eularian_data.stress(jj,ii) = max(eig([eularian_data.stressxx(jj,ii) eularian_data.stressyx(jj,ii); eularian_data.stressxy(jj,ii) eularian_data.stressyy(jj,ii)]));
            if abs(eularian_data.stress(jj,ii))> 1e8
                eularian_data.stress(jj,ii) = 0;
            end
%             Sig(logical(potentialInteractions(jj,ii,:))) = eularian_data.stress(jj,ii);
        end
    end
end
% Sig = Sig(1:N0);
% eularian_data.Sig = Sig;
% plot([Floe.poly])
warning('on',id)

end

