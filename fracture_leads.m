function Floes=fracture_leads(Floe,Nx,Ny,c2_boundary,eularian_data)
%%This function takes an input of floes and fractures each floe into a
%%specified number of smaller ones using Voronoi Tesselation
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

Floes=[]; rho_ice = 920;

x = min(c2_boundary(1,:)):(max(c2_boundary(1,:))-min(c2_boundary(1,:)))/Nx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):(max(c2_boundary(2,:))-min(c2_boundary(2,:)))/Ny:max(c2_boundary(2,:));

%Find floes and create bins
Xi=cat(1,Floe.Xi);
Yi=cat(1,Floe.Yi);
Ri=cat(1,Floe.rmax);
Binx = fix((Xi-min(x))/(max(x)-min(x))*Nx+1);
Biny = fix((Yi-min(y))/(max(y)-min(y))*Ny+1);

% Idenfity floes that are alive
SigXX = flipud(squeeze(eularian_data.stressxx)); SigYX = flipud(squeeze(eularian_data.stressyx));
SigXY = flipud(squeeze(eularian_data.stressxy)); SigYY = flipud(squeeze(eularian_data.stressyy));

for kk=1:length(Floe)
    floe = Floe(kk);
    c1=[floe.c_alpha(1,:); floe.c_alpha(2,:)];
    
    jj = Biny(kk); ii = Binx(kk);
    [V,~] = eig([SigXX(jj,ii) SigYX(jj,ii);SigXY(jj,ii) SigYY(jj,ii)]);
    r = randi([1 2],1,1);
    m = vecnorm(V(1,r)/V(2,r));
    if isinf(m)
        m = Ly;
    end
    thet = atan(m);
    r = randi([1 size(floe.interactions,1)],1,1);
    [~,r2] = min(abs(cat(1,floe.potentialInteractions.floeNum)-floe.interactions(r,1)));
    c2 = floe.potentialInteractions(r2).c-[floe.Xi; floe.Yi];
    P=InterX(c1,c2);
    count = 1;
    if ~isempty(P)
        while count < size(P,2)+1
            x1 = P(1,count); dx = 2*Ri(kk)*cos(thet);
            y1 = P(2,count); dy = 2*Ri(kk)*sin(thet);
            x2 = x1+dx; x3 = x1-dx; y2 = y1+dy; y3 = y1-dy;
            pp1=InterX(c1,[x1 x2; y1 y2]);   pp2=InterX(c1,[x1 x3; y1 y3]);
            if isempty(pp1)
                x1 = pp2(1,1); y1 = pp2(2,1); x2 = pp2(1,2); y2 = pp2(2,2);
            elseif isempty(pp2)
                x1 = pp1(1,1); y1 = pp1(2,1); x2 = pp1(1,2); y2 = pp1(2,2);
            else
                x1 = pp1(1,1); y1 = pp1(2,1); x2 = pp2(1,1); y2 = pp2(2,1);
            end
            if sqrt((x2-x1).^2-(y2-y1).^2)>10
                count = size(P,2)+1;
            end
            count = count+1;
        end
    else
        x1 = 0; x2 = 0; y1 = 0; y2 = 0;
    end
    if sqrt((x2-x1).^2-(y2-y1).^2)>10
        floe.poly = simplify(polyshape([c1(1,:);c1(2,:)]'));
        Pc1 = cutpolygon([floe.poly.Vertices;floe.poly.Vertices(1,:)],[x1, y1;x2, y2],1);
        Pc2 = cutpolygon([floe.poly.Vertices;floe.poly.Vertices(1,:)],[x1, y1;x2, y2],2);
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
            FloeNEW.angles = polyangles(FloeNEW.poly.Vertices(:,1),FloeNEW.poly.Vertices(:,2));
            FloeNEW.rmax = sqrt(max(sum((FloeNEW.poly.Vertices' - [Xi;Yi]).^2,1)));
            % n=(fix(FloeNEW.rmax/dX)+1); n=dX*(-n:n);
            % FloeNEW.Xg = n;
            % FloeNEW.Yg = n;
            % [X, Y]= meshgrid(n, n);
            % FloeNEW.X = X;
            % FloeNEW.Y = Y;
            FloeNEW.strain = floe.strain;
            FloeNEW.Stress = [0 0; 0 0];
            FloeNEW.Fx = 0; FloeNEW.Fy = 0;
            
            FloeNEW.X = FloeNEW.rmax*(2*rand(1000,1) - 1);
            FloeNEW.Y = FloeNEW.rmax*(2*rand(1000,1) - 1);
            FloeNEW.A = inpolygon(FloeNEW.X,FloeNEW.Y,FloeNEW.c_alpha(1,:),FloeNEW.c_alpha(2,:));
            % [in] = inpolygon(FloeNEW.X(:)+Xi, FloeNEW.Y(:)+Yi,FloeNEW.poly.Vertices(:,1),FloeNEW.poly.Vertices(:,2));
            % FloeNEW.A=reshape(in,length(FloeNEW.X),length(FloeNEW.X));
            
            FloeNEW.Xi = floe.Xi+Xi; FloeNEW.Yi = floe.Yi+Yi; FloeNEW.alive = 1;
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
            FloeNEW.Stress = zeros(2);
            FloeNEW.Fx = floe.Fx*area(a(p))/floe.area;
            FloeNEW.Fy = floe.Fy*area(a(p))/floe.area;
            
            Floes = [Floes FloeNEW];
            clear FloeNEW
        end
    else
        FloeNEW = floe;
        Floes = [Floes FloeNEW];
        clear FloeNEW
    end
    
end

if ~isempty(Floes)
    Floes=rmfield(Floes,{'poly'});
end

warning('on',id)
warning('on',id3)
end