function [eularian_data] = calc_eulerian_data(Floe,Nx,Ny,Nb,c2_boundary,dt,PERIODIC)
%% Function to take information of all floes and average them over a corase grained area
id = 'MATLAB:polyshape:boolOperationFailed';
warning('off',id)

%Identify only the live floes
live = cat(1,Floe.alive);
Floe(live==0)=[];
N0 = length(Floe);

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

%Create coarse grid and coarse floe variables
x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
y = fliplr(y);
dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
eularian_data.c = zeros(Ny,Nx);
eularian_data.u = zeros(Ny,Nx);
eularian_data.v = zeros(Ny,Nx);
eularian_data.du = zeros(Ny,Nx);
eularian_data.dv = zeros(Ny,Nx);
eularian_data.mom_x = zeros(Ny,Nx);
eularian_data.mom_y = zeros(Ny,Nx);
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


[xx,yy] = meshgrid(0.5*(x(1:end-1)+x(2:end)),0.5*(y(1:end-1)+y(2:end)));
xf = cat(1,Floe.Xi);
yf = cat(1,Floe.Yi);
mass = cat(1,Floe.mass);
mass(isnan(mass)==1)=0;
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

% Sig = zeros(1,length(Floe));
r_max = sqrt((dx/2)^2+(dy/2)^2);
rmax = cat(1,Floe.rmax);

%Idenfity all floes that could potentially have a piece that overlaps the
%corase areas
potentialInteractions = zeros(Ny,Nx,length(Floe));
for ii = 1:length(Floe)
    pint = sqrt((xx-xf(ii)).^2+(yy-yf(ii)).^2)-(rmax(ii)+r_max);
    pint(pint>0) = 0;
    pint(pint<0) = 1;
    potentialInteractions(:,:,ii) = pint;
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
end

if Nb > 0
    boundaries = Floe(1).poly;
    if Nb>1
        for ii = 2:Nb
            boundaries = union(boundaries,Floe(ii).poly);
        end
    end
else
    boundaries = [];
end

%% Loop to find coarse averages

for ii = 1:Nx
    for jj = 1:Ny
        if sum(logical(potentialInteractions(jj,ii,:)))>0
            %Create polyshape for the coarse area of interest
            bound = [x(ii) x(ii) x(ii+1) x(ii+1) x(ii);y(jj) y(jj+1) y(jj+1) y(jj) y(jj)];
            box = polyshape(bound(1,:), bound(2,:));
            if ~isempty(boundaries)
                box = subtract(box,boundaries);
            end
            
            %Find all floes from the potentially interacting ones that have
            %a piece in this area
            overlap = intersect(box,[Floe(logical(potentialInteractions(jj,ii,:))).poly]);
            Aover = area(overlap);
            
            
            %Calculate the concentration of floe in this area
            if sum(Aover) > 0
                polyu = union([overlap(Aover>0)]);
                Area = area(polyu);
                [~,~,zg]=meshgrid(ones(2,1),ones(2,1),Aover/Area);
                eularian_data.c(jj,ii) = Area/area(box);
            end     
        end
        
        %Now find all other variables of interest
        if eularian_data.c(jj,ii) > 0
            eularian_data.u(jj,ii) = sum(U(logical(potentialInteractions(jj,ii,:)))'.*Aover)./Area;
            eularian_data.v(jj,ii) = sum(V(logical(potentialInteractions(jj,ii,:)))'.*Aover)./Area;
            eularian_data.du(jj,ii) = sum(dU(logical(potentialInteractions(jj,ii,:)))'.*Aover)./Area;
            eularian_data.dv(jj,ii) = sum(dV(logical(potentialInteractions(jj,ii,:)))'.*Aover)./Area;
            eularian_data.mom_x(jj,ii) = sum(mass(logical(potentialInteractions(jj,ii,:)))'.*U(logical(potentialInteractions(jj,ii,:)))'.*Aover)./Area;
            eularian_data.mom_y(jj,ii) = sum(mass(logical(potentialInteractions(jj,ii,:)))'.*V(logical(potentialInteractions(jj,ii,:)))'.*Aover)./Area;
            eularian_data.force_x(jj,ii) = sum(ForceX(logical(potentialInteractions(jj,ii,:)))'.*Aover)./Area;%sum(mass(logical(potentialInteractions(jj,ii,:)))'.*dU(logical(potentialInteractions(jj,ii,:)))'.*Aover)./(dt*Area);
            eularian_data.force_y(jj,ii) = sum(ForceY(logical(potentialInteractions(jj,ii,:)))'.*Aover)./Area;%sum(mass(logical(potentialInteractions(jj,ii,:)))'.*dV(logical(potentialInteractions(jj,ii,:)))'.*Aover)./(dt*Area);
%             eularian_data.stress(jj,ii) = 0.5*trace(sum(Stress(:,:,logical(potentialInteractions(jj,ii,:))).*zg,3)./sum(Aover))/2;
            eularian_data.stressxx(jj,ii) = sum(squeeze(Stress(1,1,logical(potentialInteractions(jj,ii,:)))).*Aover')./(sum(Aover)*Area)/2;
            eularian_data.stressyx(jj,ii) = sum(squeeze(Stress(1,2,logical(potentialInteractions(jj,ii,:)))).*Aover')./(sum(Aover)*Area)/2;
            eularian_data.stressxy(jj,ii) = sum(squeeze(Stress(2,1,logical(potentialInteractions(jj,ii,:)))).*Aover')./(sum(Aover)*Area)/2;
            eularian_data.stressyy(jj,ii) = sum(squeeze(Stress(2,2,logical(potentialInteractions(jj,ii,:)))).*Aover')./(sum(Aover)*Area)/2;
            eularian_data.strainux(jj,ii) = sum(squeeze(Strain(1,1,logical(potentialInteractions(jj,ii,:)))).*Aover')./(sum(Aover)*Area)/2;
            eularian_data.strainvx(jj,ii) = sum(squeeze(Strain(1,2,logical(potentialInteractions(jj,ii,:)))).*Aover')./(sum(Aover)*Area)/2;
            eularian_data.strainuy(jj,ii) = sum(squeeze(Strain(2,1,logical(potentialInteractions(jj,ii,:)))).*Aover')./(sum(Aover)*Area)/2;
            eularian_data.strainvy(jj,ii) = sum(squeeze(Strain(2,2,logical(potentialInteractions(jj,ii,:)))).*Aover')./(sum(Aover)*Area)/2;
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
Floe=rmfield(Floe,{'poly'});

end

