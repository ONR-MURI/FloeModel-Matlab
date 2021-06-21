function Floes=fractures(floe,poly,shift)
%%This function takes an input of floes and fractures each floe into a
%%specified number of smaller ones using Voronoi Tesselation
% id ='MATLAB:polyshape:repairedBySimplify';
% warning('off',id)
% id3 ='MATLAB:polyshape:boundary3Points';
% warning('off',id3)

Floes =[];
rho_ice = 920;


for i =1:length(poly)
    R1 = poly(i);
    a = R1(area(R1)>10);
    Atot = sum(area(a));
    Mtot = floe.mass*Atot/area(polyshape(floe.c_alpha'));
    
    %%Loop through all the new shapes to calculate the new properties of
    %%each
    for p=1:length(a)
        FloeNEW.poly = rmholes(a(p));
        [Xi,Yi] = centroid(FloeNEW.poly);
        FloeNEW.area = area(FloeNEW.poly);
        FloeNEW.mass = Mtot*area(a(p))/Atot;
        FloeNEW.h = FloeNEW.mass/(rho_ice*FloeNEW.area);
        %FloeNEW.inertia_moment = PolygonMoments(FloeNEW.poly.Vertices,FloeNEW.h);
        FloeNEW.c_alpha = [(FloeNEW.poly.Vertices-[Xi Yi])' [FloeNEW.poly.Vertices(1,1)-Xi; FloeNEW.poly.Vertices(1,2)-Yi]];
        FloeNEW.c0 = FloeNEW.c_alpha;
        FloeNEW.inertia_moment = PolygonMoments(FloeNEW.c0',FloeNEW.h);
        FloeNEW.angles = polyangles(FloeNEW.poly.Vertices(:,1),FloeNEW.poly.Vertices(:,2));
        FloeNEW.rmax = sqrt(max(sum((FloeNEW.poly.Vertices' - [Xi;Yi]).^2,1)));
        FloeNEW.poly = translate(a(p),[floe.Xi,floe.Yi]);
        % n=(fix(FloeNEW.rmax/dX)+1); n=dX*(-n:n);
        % FloeNEW.Xg = n;
        % FloeNEW.Yg = n;
        % [X, Y]= meshgrid(n, n);
        % FloeNEW.X = X;
        % FloeNEW.Y = Y;
        FloeNEW.strain=floe.strain;
        FloeNEW.Stress = floe.Stress*area(a(p))/Atot;
        FloeNEW.MaxShear = floe.MaxShear*area(a(p))/Atot;
        FloeNEW.Fx = floe.Fx; FloeNEW.Fy = floe.Fy;
        FloeNEW.FxOA = 0; FloeNEW.FyOA = 0; FloeNEW.torqueOA = 0;
        
        FloeNEW.X = FloeNEW.rmax*(2*rand(1000,1) - 1);
        FloeNEW.Y = FloeNEW.rmax*(2*rand(1000,1) - 1);
        FloeNEW.A = inpolygon(FloeNEW.X,FloeNEW.Y,FloeNEW.c_alpha(1,:),FloeNEW.c_alpha(2,:));
        % [in] = inpolygon(FloeNEW.X(:)+Xi, FloeNEW.Y(:)+Yi,FloeNEW.poly.Vertices(:,1),FloeNEW.poly.Vertices(:,2));
        % FloeNEW.A=reshape(in,length(FloeNEW.X),length(FloeNEW.X));
        
        FloeNEW.Xi = floe.Xi+Xi-shift(1); FloeNEW.Yi = floe.Yi+Yi-shift(2); FloeNEW.alive = 1;
        if FloeNEW.area<1e4
            FloeNEW.alive = 0;
        end
        FloeNEW.alpha_i = 0; FloeNEW.Ui = floe.Ui; FloeNEW.Vi = floe.Vi;
        FloeNEW.dXi_p = floe.dXi_p; FloeNEW.dYi_p = floe.dYi_p;
        FloeNEW.dUi_p = floe.dUi_p; FloeNEW.dVi_p = floe.dVi_p;
        FloeNEW.dalpha_i_p = 0; FloeNEW.ksi_ice = FloeNEW.area/floe.area*floe.ksi_ice;
        FloeNEW.dksi_ice_p = floe.dksi_ice_p;
        FloeNEW.interactions = [];
        FloeNEW.potentialInteractions = [];
        FloeNEW.collision_force = 0;
        %         FloeNEW.fracture_force = 0;
        FloeNEW.collision_torque = 0;
        FloeNEW.OverlapArea = 0;
                
        Floes = [Floes FloeNEW];
        clear FloeNEW
    end
    
end

if ~isempty(Floes)
    Floes=rmfield(Floes,{'poly'});
end

% warning('on',id)
% warning('on',id3)
end