function [c,vel,accel] = calc_eulerian_data(Floe,Nx,Ny,c2_boundary)

x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));
dx = abs(x(2)-x(1));
dy = abs(y(2)-y(1));
X=cat(1,Floe.Xi);
Y=cat(1,Floe.Yi);
live = cat(1,Floe.alive);
Binx = fix((X-min(x))/(max(x)-min(x))*Nx+1);
Biny = fix((Y-min(y))/(max(y)-min(y))*Ny+1);
c = zeros(Ny,Nx);
vel.u = c;
vel.v = c;
accel.du = c;
accel.dv = c;
for ii = 1:Nx
    for jj = 1:Ny
        Area = sum(cat(1,Floe(live == 1 & Binx == ii & Biny == jj).area));
        c(jj,ii) = Area/(dx*dy);
        vel.u(jj,ii) = sum(cat(1,Floe(live == 1 & Binx == ii & Biny == jj).Ui).* ...
            cat(1,Floe(live == 1 & Binx == ii & Biny == jj).area))./Area;
        vel.v(jj,ii) = sum(cat(1,Floe(live == 1 & Binx == ii & Biny == jj).Vi).* ...
            cat(1,Floe(live == 1 & Binx == ii & Biny == jj).area))./Area;
        accel.du(jj,ii) =    sum(cat(1,Floe(live == 1 & Binx == ii & Biny == jj).dUi_p).* ...
            cat(1,Floe(live == 1 & Binx == ii & Biny == jj).area))./Area; 
        accel.dv(jj,ii) =    sum(cat(1,Floe(live == 1 & Binx == ii & Biny == jj).dVi_p).* ...
            cat(1,Floe(live == 1 & Binx == ii & Biny == jj).area))./Area; 
    end
end

end

