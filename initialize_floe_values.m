function Floe = initialize_floe_values(poly1, variance)
rho_ice=920;
poly1 = translate(poly1, [(-1)^randi([0 1])*rand (-1)^randi([0 1])*rand]);
polyout = sortregions(poly1,'area','descend');
R = regions(polyout);
poly1 = R(1);
poly1 = rmholes(poly1);
%h = 0.2*rand+0.1128;
h=2+(-1)^randi([0 1])*rand/2;
[Xi,Yi] = centroid(poly1);
    Floe.area = area(poly1);
    Floe.mass = Floe.area*h*rho_ice;
    Floe.inertia_moment = PolygonMoments(poly1.Vertices,h);
    Floe.rmax = sqrt(max(sum((poly1.Vertices' - [Xi;Yi]).^2,1)));
    Floe.Xi = Xi;
    Floe.Yi = Yi;
    Floe.alpha_i = 0;
    Floe.Ui = 0;
    Floe.Vi = 0;
    Floe.ksi_ice = 0;
    Floe.alive = 1;
    Floe.dXi_p = 0;
    Floe.dYi_p = 0;
    Floe.dUi_p = 0;
    Floe.dVi_p = 0;
    Floe.dalpha_i_p = 0;
    Floe.dksi_ice_p = 0;
    Floe.collision_force = 0;
    Floe.collision_torque = 0;
    Floe.poly = poly1;
    Floe.h = h;
    Floe.E = 1e9;
    Floe.sigma_m = 400000;
    Floe.interactions = [];
    Floe.OverlapArea = 0;
    Floe.potentialInteractions = [];

%     dX=500; % resolution of the grid inside the flow
%     n=(fix(Floe.rmax/dX)+1); n=dX*(-n:n);
%     [Xg, Yg]= meshgrid(n+Floe.Xi, n+Floe.Yi);
%     
%     % in=inpolygon(Xg(:),Yg(:),floe.poly.Vertices(:,1),floe.poly.Vertices(:,2));
%     in=inpoly2([Xg(:) Yg(:)],Floe.poly.Vertices);
%     floe_mask=reshape(in,length(Xg),length(Xg));
%     Floe.thickness = Floe.h*floe_mask+rand(length(n),length(n))*variance.*floe_mask;
    
    
%     d = logical(abs(double(floe_mask)-1));
%     d = bwdist(d)-1;
%     h = sin(pi*d)./(pi*d)+Floe.h;
%     h(isinf(1./floe_mask)) = 0;
%     h(isnan(h)) = 1 + Floe.h;
%     thickness = abs(h)+rand(length(n),length(n))*(min(abs(h(floe_mask))));
%     thickness(isinf(1./floe_mask))=0;
    
end