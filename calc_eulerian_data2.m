function [c,vel,accel] = calc_eulerian_data2(Floe,Nx,Ny,c2_boundary)

x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
y = fliplr(y);
c = zeros(Ny,Nx);
vel.u = c;
vel.v = c;
accel.du = c;
accel.dv = c;
for ii = 1:Nx
    for jj = 1:Ny
        bound = [x(ii) x(ii) x(ii+1) x(ii+1) x(ii);y(jj) y(jj+1) y(jj+1) y(jj) y(jj)];
        box = polyshape(bound(1,:), bound(2,:));
        overlap = intersect(box,[Floe.poly]);
        Aover = area(overlap);
        Area = sum(Aover);
        c(jj,ii) = Area/area(box);
        if c(jj,ii) > 0
            vel.u(jj,ii) = sum(cat(1,Floe.Ui)'.*Aover)./Area;
            vel.v(jj,ii) = sum(cat(1,Floe.Vi)'.*Aover)./Area;
            accel.du(jj,ii) = sum(cat(1,Floe.dUi_p)'.*Aover)./Area;
            accel.dv(jj,ii) = sum(cat(1,Floe.dVi_p)'.*Aover)./Area;
        end
    end
end
end

