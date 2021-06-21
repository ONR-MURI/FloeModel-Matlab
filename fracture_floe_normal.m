function Floes=fracture_floe_normal(Floe,m)
%%This function takes an input of floes and fractures each floe into a
%%specified number of smaller ones using Voronoi Tesselation
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

Floes=[]; rho_ice = 920;

for kk=1:length(Floe)
if kk >1 && isempty(Floes)
xx = 1; xx(1) = [1 2];
end

    floe=Floe(kk);
%     if ~isempty(floe.interactions)
    if ~isempty(floe.potentialInteractions)
        c1 = floe.c_alpha+[floe.Xi;floe.Yi];
        contacts =[];
        for ii = 1:length(floe.potentialInteractions)
            P = InterX(c1,floe.potentialInteractions(ii).c);
            contacts = [contacts P];
        end
        if ~isempty(contacts)
            FracPoint = contacts(:,randi([1 length(contacts)]));
%         a = floe.interactions(:,4:5);
%         center = mean(a-[floe.Xi floe.Yi]);
%         FracPoint = 2*(rand(2,1)-0.5).*(floe.rmax-sign(center).*abs(center)');
%         thet2 = atan(2*floe.Stress(2)/(floe.Stress(1)-floe.Stress(4)));
%         m = [tan(thet2/2) -tan(thet2/2)];%[tan(thet2/2) tan(pi/2+thet2/2)];
            ii = randi([1 2]);
            dir = m(ii);
            if abs(dir) > 1e3
                dir = sign(m(ii))*1e3;
            end
%         floe.poly = simplify(polyshape(floe.c_alpha'+[floe.Xi floe.Yi]));
%         Pc1 = cutpolygon([floe.poly.Vertices;floe.poly.Vertices(1,:)],[FracPoint(1)+floe.Xi, FracPoint(2)+floe.Yi;FracPoint(1)+floe.Xi+1, FracPoint(2)+floe.Yi+dir],1);
%         Pc2 = cutpolygon([floe.poly.Vertices;floe.poly.Vertices(1,:)],[FracPoint(1)+floe.Xi, FracPoint(2)+floe.Yi;FracPoint(1)+floe.Xi+1, FracPoint(2)+floe.Yi+dir],2);
            floe.poly = simplify(polyshape(c1'));
            Pc1 = cutpolygon([floe.poly.Vertices;floe.poly.Vertices(1,:)],[FracPoint(1), FracPoint(2);FracPoint(1)+1, FracPoint(2)+dir],1);
            Pc2 = cutpolygon([floe.poly.Vertices;floe.poly.Vertices(1,:)],[FracPoint(1), FracPoint(2);FracPoint(1)+1, FracPoint(2)+dir],2);
            poly1 = polyshape(Pc1(1:end-1,:));
            poly2 = polyshape(Pc2(1:end-1,:));
            R1 = [regions(poly1); regions(poly2)];
            R1 = translate(R1,-[floe.Xi floe.Yi]);
            a = R1(area(R1)>1e3);
        else
            a = 0;
            floe.poly = simplify(polyshape(floe.c_alpha'+[floe.Xi floe.Yi]));
        end
    else
        a = 0;
        floe.poly = simplify(polyshape(floe.c_alpha'+[floe.Xi floe.Yi]));
    end
    
    %%Loop through all the new shapes to calculate the new properties of
    %%each
    if length(a) > 1
        for p=1:length(a)
            FloeNEW.poly = rmholes(a(p));
            [Xi,Yi] = centroid(FloeNEW.poly);
            FloeNEW.area = area(FloeNEW.poly);
            FloeNEW.mass = floe.mass*area(a(p))/floe.area;
            FloeNEW.h = floe.mass*area(a(p))/(rho_ice*FloeNEW.area*floe.area);
            FloeNEW.c_alpha = [(FloeNEW.poly.Vertices-[Xi Yi])' [FloeNEW.poly.Vertices(1,1)-Xi; FloeNEW.poly.Vertices(1,2)-Yi]];
            FloeNEW.c0 = FloeNEW.c_alpha;
            FloeNEW.inertia_moment = PolygonMoments(FloeNEW.c0',FloeNEW.h);
            FloeNEW.angles = polyangles(FloeNEW.poly.Vertices(:,1),FloeNEW.poly.Vertices(:,2));
            FloeNEW.rmax = sqrt(max(sum((FloeNEW.poly.Vertices' - [Xi;Yi]).^2,1)));
            FloeNEW.strain = floe.strain;
            FloeNEW.Stress = floe.Stress*area(a(p))/floe.area;
            FloeNEW.MaxShear = floe.MaxShear*area(a(p))/floe.area;
            FloeNEW.Fx = floe.Fx; FloeNEW.Fy = floe.Fy;
            FloeNEW.FxOA = 0; FloeNEW.FyOA = 0; FloeNEW.torqueOA = 0;
            
            err = 1;
            while err > 0.1
                FloeNEW.X = FloeNEW.rmax*(2*rand(1000,1) - 1);
                FloeNEW.Y = FloeNEW.rmax*(2*rand(1000,1) - 1);
                FloeNEW.A = inpolygon(FloeNEW.X,FloeNEW.Y,FloeNEW.c_alpha(1,:),FloeNEW.c_alpha(2,:));
                err = (sum(FloeNEW.A)/1000*4*FloeNEW.rmax^2-FloeNEW.area)/FloeNEW.area;
            end
            
            FloeNEW.Xi = floe.Xi+Xi; FloeNEW.Yi = floe.Yi+Yi; FloeNEW.alive = 1;
            FloeNEW.alpha_i = 0; FloeNEW.Ui = floe.Ui; FloeNEW.Vi = floe.Vi;
            FloeNEW.dXi_p = floe.dXi_p; FloeNEW.dYi_p = floe.dYi_p;
            FloeNEW.dUi_p = floe.dUi_p; FloeNEW.dVi_p = floe.dVi_p;
            FloeNEW.dalpha_i_p = 0; FloeNEW.ksi_ice = floe.ksi_ice;%FloeNEW.area/floe.area*floe.ksi_ice;
            FloeNEW.dksi_ice_p = floe.dksi_ice_p;
            FloeNEW.interactions = [];%floe.interactions;
            FloeNEW.potentialInteractions = [];%floe.potentialInteractions;
            FloeNEW.collision_force = 0;
            FloeNEW.collision_torque = 0;
            FloeNEW.OverlapArea = 0;
            FloeNEW.Stress = zeros(2);
            FloeNEW.Fx = floe.Fx*area(a(p))/floe.area;
            FloeNEW.Fy = floe.Fy*area(a(p))/floe.area;
            
            Floes = [Floes FloeNEW];
        end
        clear FloeNEW
    else
        Floes = [Floes floe];
        clear FloeNEW;
    end
end
    
if ~isempty(Floes)
    Floes=rmfield(Floes,{'poly'});
end

warning('on',id)
warning('on',id3)
end