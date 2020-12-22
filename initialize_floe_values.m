function FloeNEW = initialize_floe_values(poly1, height)
%%This function populates all the fields of the floes upon their creation

rho_ice=920;
dX = 480;

polyout = sortregions(poly1,'area','descend');
R = regions(polyout);
poly1 = R(1);
polya = rmholes(poly1);
h=height.mean+(-1)^randi([0 1])*rand*height.delta;

FloeNEW.poly = poly1;
[Xi,Yi] = centroid(FloeNEW.poly);
FloeNEW.area = area(FloeNEW.poly);
FloeNEW.h = h;
FloeNEW.mass = FloeNEW.area*h*rho_ice;
FloeNEW.inertia_moment = PolygonMoments(polya.Vertices,h);
FloeNEW.c_alpha = [(poly1.Vertices-[Xi Yi])' [poly1.Vertices(1,1)-Xi; poly1.Vertices(1,2)-Yi]];
FloeNEW.c0 = FloeNEW.c_alpha;
FloeNEW.angles = polyangles(polya.Vertices(:,1),polya.Vertices(:,2));
FloeNEW.rmax = sqrt(max(sum((FloeNEW.poly.Vertices' - [Xi;Yi]).^2,1)));
% n=(fix(FloeNEW.rmax/dX)+1); n=dX*(-n:n);
% FloeNEW.Xg = n;
% FloeNEW.Yg = n;
% [X, Y]= meshgrid(n, n);
% FloeNEW.X = X;
% FloeNEW.Y = Y;
FloeNEW.Stress = [0 0; 0 0];
FloeNEW.strain = [0 0; 0 0];
FloeNEW.Fx = 0; FloeNEW.Fy = 0;

FloeNEW.X = FloeNEW.rmax*(2*rand(1000,1) - 1);
FloeNEW.Y = FloeNEW.rmax*(2*rand(1000,1) - 1);
FloeNEW.A = inpolygon(FloeNEW.X,FloeNEW.Y,FloeNEW.c_alpha(1,:),FloeNEW.c_alpha(2,:));
% [in] = inpolygon(FloeNEW.X(:)+Xi, FloeNEW.Y(:)+Yi,FloeNEW.poly.Vertices(:,1),FloeNEW.poly.Vertices(:,2));
% FloeNEW.A=reshape(in,length(FloeNEW.X),length(FloeNEW.X));

FloeNEW.Xi = Xi; FloeNEW.Yi = Yi; FloeNEW.alive = 1;
FloeNEW.alpha_i = 0; FloeNEW.Ui = 0; FloeNEW.Vi = 0;
FloeNEW.dXi_p = 0; FloeNEW.dYi_p = 0;
FloeNEW.dUi_p = 0; FloeNEW.dVi_p = 0;
FloeNEW.dalpha_i_p = 0; FloeNEW.ksi_ice = 0;
FloeNEW.dksi_ice_p = 0;
FloeNEW.interactions = [];
FloeNEW.potentialInteractions = [];
FloeNEW.collision_force = 0;
%         FloeNEW.fracture_force = 0;
FloeNEW.collision_torque = 0;
FloeNEW.OverlapArea = 0;

if isnan(FloeNEW.inertia_moment)
    xx = 1;
    xx(1) = [1 2];
end

end