
function Floes=fracture_floe(Floe,N)

Floes=[]; k=1; rho_ice = 920;

for kk=1:length(Floe)
    %% 
    
    floe=Floe(kk);


X = floe.Xi+floe.rmax*(2*rand(N,1)-1);
Y = floe.Yi+floe.rmax*(2*rand(N,1)-1);


boundingbox=[-1 ,-1; 1,-1; 1,1; -1 ,1]*floe.rmax+[floe.Xi floe.Yi];
[~, b] = polybnd_voronoi([X Y],boundingbox);

%% 

for i=1:length(b)
    a=regions(intersect(floe.poly,polyshape(b{i}))); 
    for p=1:length(a)
        if k==1 
            Floes.poly = a(p);
            [Xi,Yi] = centroid(Floes.poly);
            Floes.area = area(Floes.poly);
            Floes.rmax = sqrt(max(sum((Floes.poly.Vertices' - [Xi;Yi]).^2,1)));
            Floes.Xi = Xi; Floes.Yi = Yi; Floes.alive = 1;
            Floes.collision_force = 0;
            Floes.collision_torque = 0;
            Floes.h = floe.h;
            Floes.E = 1e9;
            Floes.sigma_m = 400000;
            Floes.interactions = [];
            Floes.OverlapArea = 0;
            Floes.potentialInteractions = [];
            Floes.alpha_i = floe.alpha_i; Floes.Ui = floe.Ui; Floes.Vi = floe.Vi;
            Floes.dXi_p = floe.dXi_p; Floes.dYi_p = floe.dYi_p;
            Floes.dUi_p = floe.dUi_p; Floes.dVi_p = floe.dVi_p;
            Floes.dalpha_i_p = floe.dalpha_i_p; Floes.ksi_ice = Floes.area/floe.area*floe.ksi_ice;
            Floes.dksi_ice_p = Floes.area/floe.area*floe.dksi_ice_p;
            Floes.vorX = floe.vorX; Floes.vorY = floe.vorY; Floes.vorbox = floe.vorbox;
            j = 1; 
            for jj = 1:length(floe.SubFloes)
                polyout = intersect(floe.SubFloes(jj).poly,Floes.poly);
                if area(polyout) > 0
                    R = regions(polyout);
                    polynew = R(1);
                    polyout = rmholes(polynew);
                    Floes.SubFloes(j).poly = polyout;
                    Floes.SubFloes(j).h = floe.SubFloes(jj).h;
                    areaS(j) = area(Floes.SubFloes(j).poly);
                    inertia(j) = PolygonMoments(Floes.SubFloes(j).poly.Vertices,Floes.SubFloes(j).h);
                    [Xi,Yi] = centroid(Floes.SubFloes(j).poly);
                    centers(j,:) = [Xi,Yi];
                    j = j+1;
                end
            end
            if j > 1
                Floes.mass = sum(rho_ice*areaS'.*cat(1,Floes.SubFloes.h));
                Floes.Xm = sum(rho_ice*areaS'.*cat(1,Floes.SubFloes.h).*centers(:,1))./Floes.mass;
                Floes.Ym = sum(rho_ice*areaS'.*cat(1,Floes.SubFloes.h).*centers(:,2))./Floes.mass;
                Floes.inertia_moment = sum(inertia'+cat(1,Floes.SubFloes.h).*sqrt((centers(:,1)-Floes.Xm).^2+(centers(:,2)-Floes.Ym).^2));
                k=k+1;
            else
                Floes = [];
            end
        else
            FloeNEW.poly = a(p);
            [Xi,Yi] = centroid(FloeNEW.poly);
            FloeNEW.area = area(FloeNEW.poly);
            FloeNEW.rmax = sqrt(max(sum((FloeNEW.poly.Vertices' - [Xi;Yi]).^2,1)));
            FloeNEW.Xi = Xi; FloeNEW.Yi = Yi; FloeNEW.alive = 1;
            FloeNEW.collision_force = 0;
            FloeNEW.collision_torque = 0;
            FloeNEW.h = floe.h;
            FloeNEW.E = 1e9;
            FloeNEW.sigma_m = 400000;
            FloeNEW.interactions = [];
            FloeNEW.OverlapArea = 0;
            FloeNEW.potentialInteractions = [];
            FloeNEW.alpha_i = floe.alpha_i; FloeNEW.Ui = floe.Ui; FloeNEW.Vi = floe.Vi;
            FloeNEW.dXi_p = floe.dXi_p; FloeNEW.dYi_p = floe.dYi_p;
            FloeNEW.dUi_p = floe.dUi_p; FloeNEW.dVi_p = floe.dVi_p;
            FloeNEW.dalpha_i_p = floe.dalpha_i_p; FloeNEW.ksi_ice = FloeNEW.area/floe.area*floe.ksi_ice;
            FloeNEW.dksi_ice_p = FloeNEW.area/floe.area*floe.dksi_ice_p;
            FloeNEW.vorX = floe.vorX; FloeNEW.vorY = floe.vorY; FloeNEW.vorbox = floe.vorbox;
            j = 1;
            clear areaS; clear inertia; clear centers;
            for jj = 1:length(floe.SubFloes)
                polyout = intersect(floe.SubFloes(jj).poly,FloeNEW.poly);
                if area(polyout) > 0
                    R = regions(polyout);
                    polynew = R(1);
                    polyout = rmholes(polynew);
                    FloeNEW.SubFloes(j).poly = polyout;
                    FloeNEW.SubFloes(j).h = floe.SubFloes(jj).h;
                    areaS(j) = area(FloeNEW.SubFloes(j).poly);
                    inertia(j) = PolygonMoments(FloeNEW.SubFloes(j).poly.Vertices,FloeNEW.SubFloes(j).h);
                    [Xi,Yi] = centroid(FloeNEW.SubFloes(j).poly);
                    centers(j,:) = [Xi,Yi];
                    j = j+1;
                end
            end
            if j > 1
                FloeNEW.mass = sum(rho_ice*areaS'.*cat(1,FloeNEW.SubFloes.h));
                FloeNEW.Xm = sum(rho_ice*areaS'.*cat(1,FloeNEW.SubFloes.h).*centers(:,1))./FloeNEW.mass;
                FloeNEW.Ym = sum(rho_ice*areaS'.*cat(1,FloeNEW.SubFloes.h).*centers(:,2))./FloeNEW.mass;
                FloeNEW.inertia_moment = sum(inertia'+cat(1,FloeNEW.SubFloes.h).*sqrt((centers(:,1)-FloeNEW.Xm).^2+(centers(:,2)-FloeNEW.Ym).^2));
                Floes(k) = FloeNEW;
                k=k+1;
            end
            clear FloeNEW
        end

    end
    
end




end

% figure; 
% plot(floe.poly)
% hold on;
% plot(X,Y,'x');
% plot([Floes.poly])



end