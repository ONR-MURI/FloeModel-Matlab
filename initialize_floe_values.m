function Floe = initialize_floe_values(poly1, SUBFLOES,PACKING)

rho_ice=920;
poly1 = translate(poly1, [(-1)^randi([0 1])*rand (-1)^randi([0 1])*rand]);
polyout = sortregions(poly1,'area','descend');
R = regions(polyout);
poly1 = R(1);
poly1 = rmholes(poly1);
%h = 0.2*rand+0.1128;
h=2;%+(-1)^randi([0 1])*rand/2;
if PACKING
    h = 0.15;
end
[Xi,Yi] = centroid(poly1);
    Floe.area = area(poly1);
    Floe.mass = 0;%Floe.area*h*rho_ice;
    Floe.inertia_moment = 0;%PolygonMoments(poly1.Vertices,h);
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
    
    
    EXISTING = false;
    N1 = ceil(sqrt(Floe.area/1e6))+2;
    [subFloes,Floe.vorX,Floe.vorY,Floe.vorbox] = create_subfloes(Floe,N1,EXISTING);
    if ~SUBFLOES
        subFloes = Floe.poly;
    end
    N = length(subFloes);
    areaS = zeros(N,1);
    inertia = zeros(N,1);
    centers = zeros(N,2);
    for ii = 1:N
        R = regions(subFloes(ii));
        polynew = R(1);
        Floe.SubFloes(ii).poly = rmholes(polynew);
        Floe.SubFloes(ii).h = h;
        areaS(ii) = area(Floe.SubFloes(ii).poly);
        if areaS(ii) == 0
            x = 1;
            x(1) = [1 2];
        end
        inertia(ii) = PolygonMoments(Floe.SubFloes(ii).poly.Vertices,Floe.SubFloes(ii).h);
        [Xi,Yi] = centroid(Floe.SubFloes(ii).poly);
        centers(ii,:) = [Xi,Yi];
    end
    Floe.mass = sum(rho_ice*areaS.*cat(1,Floe.SubFloes.h));
    Floe.Xm = sum(rho_ice*areaS.*cat(1,Floe.SubFloes.h).*centers(:,1))./Floe.mass;
    Floe.Ym = sum(rho_ice*areaS.*cat(1,Floe.SubFloes.h).*centers(:,2))./Floe.mass;
    Floe.inertia_moment = sum(inertia+cat(1,Floe.SubFloes.h).*sqrt((centers(:,1)-Floe.Xm).^2+(centers(:,2)-Floe.Ym).^2));
    Floe.rmax = sqrt(max(sum((Floe.poly.Vertices' - [Floe.Xi; Floe.Yi]).^2,1)));
    
 
end