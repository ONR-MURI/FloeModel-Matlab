function Floe = initialize_floe_values(poly1, height,SHIFT, SUBFLOES)
%%This function populates all the fields of the floes upon their creation

rho_ice=920;
if SHIFT
    poly1 = translate(poly1, [(-1)^randi([0 1])*rand (-1)^randi([0 1])*rand]);
end
polyout = sortregions(poly1,'area','descend');
R = regions(polyout);
poly1 = R(1);
h=height.mean+(-1)^randi([0 1])*rand*height.delta;
[Xi,Yi] = centroid(poly1);
Floe.area = area(poly1);
Floe.mass = 0;%This will be filled in later
Floe.inertia_moment = 0;%This will be filled in later
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

%Now account for potential variable thickness for different regions of the
%flow and calculate all properties that this variation might impact here by
%looping through reach regoin
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
    if Floe.SubFloes(ii).poly.NumHoles > 0
        breaks = isnan(polynew.Vertices(:,1));
        I = find(breaks == 1);
        I = [0 I' length(breaks)+1];
        inertia(ii) = 0;
        for jj = length(I) -1
            inertia(ii) = inertia(ii) + PolygonMoments(Floe.SubFloes(ii).poly.Vertices(I(jj)+1:I(jj+1)-1,:),Floe.SubFloes(ii).h);
        end
    else
        inertia(ii) = PolygonMoments(Floe.SubFloes(ii).poly.Vertices,Floe.SubFloes(ii).h);
    end
    [Xi,Yi] = centroid(Floe.SubFloes(ii).poly);
    centers(ii,:) = [Xi,Yi];
end
Floe.mass = sum(rho_ice*areaS.*cat(1,Floe.SubFloes.h));
Floe.Xm = sum(rho_ice*areaS.*cat(1,Floe.SubFloes.h).*centers(:,1))./Floe.mass;
Floe.Ym = sum(rho_ice*areaS.*cat(1,Floe.SubFloes.h).*centers(:,2))./Floe.mass;
Floe.inertia_moment = sum(inertia+cat(1,Floe.SubFloes.h).*sqrt((centers(:,1)-Floe.Xm).^2+(centers(:,2)-Floe.Ym).^2));
Floe.rmax = sqrt(max(sum((Floe.poly.Vertices' - [Floe.Xi; Floe.Yi]).^2,1)));

if isnan(Floe.inertia_moment)
    xx = 1;
    xx(1) = [1 2];
end
end