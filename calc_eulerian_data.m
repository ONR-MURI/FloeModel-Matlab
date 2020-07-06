function [eularian_data] = calc_eulerian_data(Floe,Nx,Ny,c2_boundary,PERIODIC)
%% Function to take information of all floes and average them over a corase grained area
id = 'MATLAB:polyshape:boolOperationFailed';
warning('off',id)

%Identify only the live floes
live = cat(1,Floe.alive);
Floe(live==0)=[];

%Create ghost floes for periodic floe states
Lx= max(c2_boundary(1,:));
Ly= max(c2_boundary(2,:));%c2 must be symmetric around x=0 for channel boundary conditions.
if PERIODIC
    
    ghostFloeX=[];
    ghostFloeY=[];
    parent=[];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(Floe(i).poly.Vertices(:,1)))>Lx)
            
            ghostFloeX=[ghostFloeX  Floe(i)];
            ghostFloeX(end).poly=translate(Floe(i).poly,[-2*Lx*sign(x(i)) 0]);
            for ii = 1:length(Floe(i).SubFloes)
                ghostFloeX(end).SubFloes(ii).poly=translate(Floe(i).SubFloes(ii).poly,[-2*Lx*sign(x(i)) 0]);
            end
            ghostFloeX(end).Xi=Floe(i).Xi-2*Lx*sign(x(i));
            ghostFloeX(end).vorX=Floe(i).vorX-2*Lx*sign(x(i));
            ghostFloeX(end).vorbox(:,1)=Floe(i).vorbox(:,1)-2*Lx*sign(x(i));
            parent=[parent  i];
            
        end
        
        
    end
    
    Floe=[Floe ghostFloeX];
    
    x=cat(1,Floe.Xi);
    y=cat(1,Floe.Yi);
    alive=cat(1,Floe.alive);
    
    for i=1:length(Floe)
        
        %   if alive(i) && (x(i)>Lx-rmax(i)) || (x(i)<-Lx+rmax(i))
        if alive(i) && (max(abs(Floe(i).poly.Vertices(:,2)))>Ly)
            
            ghostFloeY=[ghostFloeY  Floe(i)];
            ghostFloeY(end).poly=translate(Floe(i).poly,[0 -2*Ly*sign(y(i))]);
            for ii = 1:length(Floe(i).SubFloes)
                ghostFloeY(end).SubFloes(ii).poly=translate(Floe(i).SubFloes(ii).poly,[0 -2*Ly*sign(y(i))]);
            end
            ghostFloeY(end).Yi=Floe(i).Yi-2*Ly*sign(y(i));
            ghostFloeY(end).vorY=Floe(i).vorY-2*Ly*sign(y(i));
            ghostFloeY(end).vorbox(:,2)=Floe(i).vorbox(:,2)-2*Ly*sign(y(i));
            parent=[parent  i];
            
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
[xx,yy] = meshgrid(0.5*(x(1:end-1)+x(2:end)),0.5*(y(1:end-1)+y(2:end)));
xf = cat(1,Floe.Xi);
yf = cat(1,Floe.Yi);
mass = cat(1,Floe.mass);
mass(isnan(mass)==1)=0;
U = cat(1,Floe.Ui);
U(isnan(U)==1)=0;
V = cat(1,Floe.Vi);
V(isnan(V)==1)=0;
dU = cat(1,Floe.dUi_p);
dU(isnan(dU)==1)=0;
dV = cat(1,Floe.dVi_p);
dV(isnan(dV)==1)=0;
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
end
%% Loop to find coarse averages

for ii = 1:Nx
    for jj = 1:Ny
        if sum(logical(potentialInteractions(jj,ii,:)))>0
            %Create polyshape for the coarse area of interest
            bound = [x(ii) x(ii) x(ii+1) x(ii+1) x(ii);y(jj) y(jj+1) y(jj+1) y(jj) y(jj)];
            box = polyshape(bound(1,:), bound(2,:));
            
            %Find all floes from the potentially interacting ones that have
            %a piece in this area
            overlap = intersect(box,[Floe(logical(potentialInteractions(jj,ii,:))).poly]);
            Aover = area(overlap);
            
            %Calculate the concentration of floe in this area
            if sum(Aover) > 0
                polyu = union([overlap(Aover>0)]);
                Area = area(polyu);
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
            eularian_data.force_x(jj,ii) = sum(mass(logical(potentialInteractions(jj,ii,:)))'.*dU(logical(potentialInteractions(jj,ii,:)))'.*Aover)./Area;
            eularian_data.force_y(jj,ii) = sum(mass(logical(potentialInteractions(jj,ii,:)))'.*dV(logical(potentialInteractions(jj,ii,:)))'.*Aover)./Area;
        end
    end
end

warning('on',id)
end

