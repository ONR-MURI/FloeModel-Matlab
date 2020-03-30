function [eularian_data] = calc_eulerian_data2(Floe,Nx,Ny,c2_boundary)

x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
y = fliplr(y);
dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
eularian_data.c = zeros(Ny,Nx);
eularian_data.vel.u = eularian_data.c;
eularian_data.vel.v = eularian_data.c;
eularian_data.accel.du = eularian_data.c;
eularian_data.accel.dv = eularian_data.c;
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
potentialInteractions = zeros(Ny,Nx,length(Floe));
for ii = 1:length(Floe)
    pint = sqrt((xx-xf(ii)).^2+(yy-yf(ii)).^2)-(rmax(ii)+r_max);
    pint(pint>0) = 0;
    pint(pint<0) = 1;
    potentialInteractions(:,:,ii) = pint;
end
for ii = 1:Nx
    for jj = 1:Ny
        if sum(logical(potentialInteractions(jj,ii,:)))>0
            bound = [x(ii) x(ii) x(ii+1) x(ii+1) x(ii);y(jj) y(jj+1) y(jj+1) y(jj) y(jj)];
            box = polyshape(bound(1,:), bound(2,:));
            overlap = intersect(box,[Floe(logical(potentialInteractions(jj,ii,:))).poly]);
            Aover = area(overlap);
            Area = sum(Aover);
            eularian_data.c(jj,ii) = Area/area(box);
        end
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
end

