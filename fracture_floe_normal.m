function Floes=fracture_floe_normal(Floe,cf)
%%This function takes an input of floes and fractures each floe into a
%%specified number of smaller ones using Voronoi Tesselation
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

Floes=[]; rho_ice = 920;
normal = cf/vecnorm(cf);

for kk=1:length(Floe)
    floe=Floe(kk);
    if abs(normal(1)) < 1e-3
        cf = [ 0.9 0.01];
        dir = cf/vecnorm(cf);
        normal = [-dir(2) dir(1)];
    end
    floe.poly = simplify(polyshape(floe.c_alpha'+[floe.Xi floe.Yi]));
    Pc1 = cutpolygon([floe.poly.Vertices;floe.poly.Vertices(1,:)],[floe.Xi, floe.Yi;floe.Xi+normal(1), floe.Yi+normal(2)],1);
    Pc2 = cutpolygon([floe.poly.Vertices;floe.poly.Vertices(1,:)],[floe.Xi, floe.Yi;floe.Xi+normal(1), floe.Yi+normal(2)],2);
    poly1 = polyshape(Pc1(1:end-1,:));
    poly2 = polyshape(Pc2(1:end-1,:));
    R1 = [regions(poly1); regions(poly2)];
    a = R1(area(R1)>1e3);
    
    %%Loop through all the new shapes to calculate the new properties of
    %%each
    for p=1:length(a)
        FloeNEW.poly = rmholes(a(p));
        [Xi,Yi] = centroid(FloeNEW.poly);
        FloeNEW.area = area(FloeNEW.poly);
        FloeNEW.mass = floe.mass*area(a(p))/floe.area;
        FloeNEW.h = floe.mass*area(a(p))/(rho_ice*FloeNEW.area*floe.area);
        FloeNEW.inertia_moment = PolygonMoments(FloeNEW.poly.Vertices,FloeNEW.h);
        FloeNEW.c_alpha = [(FloeNEW.poly.Vertices-[Xi Yi])' [FloeNEW.poly.Vertices(1,1)-Xi; FloeNEW.poly.Vertices(1,2)-Yi]];
        FloeNEW.c0 = FloeNEW.c_alpha;
        FloeNEW.rmax = sqrt(max(sum((FloeNEW.poly.Vertices' - [Xi;Yi]).^2,1)));
        FloeNEW.Xg = floe.Xg;
        FloeNEW.Yg = floe.Yg;
        FloeNEW.X = floe.X;
        FloeNEW. Y = floe.Y;
        
        [in] = inpolygon(FloeNEW.X(:), FloeNEW.Y(:),FloeNEW.poly.Vertices(:,1),FloeNEW.poly.Vertices(:,2));
        floe_mask=reshape(in,length(FloeNEW.X),length(FloeNEW.X));
        FloeNEW.A = flipud(floe_mask);
        
        FloeNEW.Xi = Xi; FloeNEW.Yi = Yi; FloeNEW.alive = 1;
        FloeNEW.alpha_i = 0; FloeNEW.Ui = floe.Ui; FloeNEW.Vi = floe.Vi;
        FloeNEW.dXi_p = floe.dXi_p; FloeNEW.dYi_p = floe.dYi_p;
        FloeNEW.dUi_p = floe.dUi_p; FloeNEW.dVi_p = floe.dVi_p;
        FloeNEW.dalpha_i_p = 0; FloeNEW.ksi_ice = FloeNEW.area/floe.area*floe.ksi_ice;
        FloeNEW.dksi_ice_p = FloeNEW.area/floe.area*floe.dksi_ice_p;
        FloeNEW.interactions = [];
        FloeNEW.potentialInteractions = [];
        FloeNEW.collision_force = 0;
%         FloeNEW.fracture_force = 0;
        FloeNEW.collision_torque = 0;
%         FloeNEW.OverlapArea = 0;
%         FloeNEW.Stress = zeros(2);
        
        Floes = [Floes FloeNEW];
    end
    clear FloeNEW
end

if ~isempty(Floes)
    Floes=rmfield(Floes,{'poly'});
end

warning('on',id)
warning('on',id3)
end