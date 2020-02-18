function [c,vel,accel] = calc_eulerian_data2(Floe,Nx,Ny,c2_boundary)

x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
y = fliplr(y);
dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
c = zeros(Ny,Nx);
vel.u = c;
vel.v = c;
accel.du = c;
accel.dv = c;
[xx,yy] = meshgrid(0.5*(x(1:end-1)+x(2:end)),0.5*(y(1:end-1)+y(2:end)));
xf = cat(1,Floe.Xi);
yf = cat(1,Floe.Yi);
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
        bound = [x(ii) x(ii) x(ii+1) x(ii+1) x(ii);y(jj) y(jj+1) y(jj+1) y(jj) y(jj)];
        box = polyshape(bound(1,:), bound(2,:));
        overlap = intersect(box,[Floe(logical(potentialInteractions(jj,ii,:))).poly]);
        Aover = area(overlap);
        Area = sum(Aover);
        c(jj,ii) = Area/area(box);
        if c(jj,ii) > 0
            vel.u(jj,ii) = sum(cat(1,Floe(logical(potentialInteractions(jj,ii,:))).Ui)'.*Aover)./Area;
            vel.v(jj,ii) = sum(cat(1,Floe(logical(potentialInteractions(jj,ii,:))).Vi)'.*Aover)./Area;
            accel.du(jj,ii) = sum(cat(1,Floe(logical(potentialInteractions(jj,ii,:))).dUi_p)'.*Aover)./Area;
            accel.dv(jj,ii) = sum(cat(1,Floe(logical(potentialInteractions(jj,ii,:))).dVi_p)'.*Aover)./Area;
        end
    end
end
end

