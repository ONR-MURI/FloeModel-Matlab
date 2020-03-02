function [Vd] = calc_vol_dissolved(Floe,Nx,Ny,c2_boundary_poly)

x = min(c2_boundary_poly.Vertices(1,:)):(max(c2_boundary_poly.Vertices(1,:))-min(c2_boundary_poly.Vertices(1,:)))/Nx:max(c2_boundary_poly.Vertices(1,:));
y = min(c2_boundary_poly.Vertices(2,:)):(max(c2_boundary_poly.Vertices(2,:))-min(c2_boundary_poly.Vertices(2,:)))/Ny:max(c2_boundary_poly.Vertices(2,:));
X=cat(1,Floe.Xi);
Y=cat(1,Floe.Yi);
live = cat(1,Floe.alive);
Binx = fix((X-min(x))/(max(x)-min(x))*Nx+1);
Biny = fix((Y-min(y))/(max(y)-min(y))*Ny+1);
Vd = zeros(Ny,Nx);
for ii = 1:Nx
    for jj = 1:Ny
        Vd(jj,ii) = sum(cat(1,Floe(live == 1 & Binx == ii & Biny == jj).mass));%.* ...
            %cat(1,Floe(live == 1 & Binx == ii & Biny == jj).area).*...
            %cat(1,Floe(live == 1 & Binx == ii & Biny == jj).h));
    end
end
Vd = flipud(Vd);
end

