function [FloeWhole,FloeNEW] = ridge_values_update(FloeWhole,FloeLess, overlap1, overlap2, V, Area)
%%This function calculates the values for the floes after their shapes and
%%mass get updated after ridging

%Store a copy of this structure for later reference
Floe2 = FloeLess;

%   Calculate values for Floe gaining the new ice
rho_ice = 920;
Nw = length(FloeWhole.SubFloes);
for ii = 1:length(overlap1)
    FloeWhole.SubFloes(overlap1(ii)).h = FloeWhole.SubFloes(overlap1(ii)).h + V/Area;
    if FloeWhole.SubFloes(overlap1(ii)).h > 30
        FloeWhole.SubFloes(overlap1(ii)).h = 30;
    end
end
FloeWhole.mass = FloeWhole.mass+V*rho_ice;
areaS = zeros(Nw,1);
inertia = zeros(Nw,1);
centers = zeros(Nw,2);
for ii = 1:Nw
    areaS(ii) = area(FloeWhole.SubFloes(ii).poly);
    if FloeWhole.SubFloes(ii).poly.NumHoles > 0
        breaks = isnan(FloeWhole.SubFloes(ii).poly.Vertices(:,1));
        I = find(breaks == 1);
        I = [0 I' length(breaks)+1];
        inertia(ii) = 0;
        for kk = length(I) -1
            inertia(ii) = inertia(ii) + PolygonMoments(FloeWhole.SubFloes(ii).poly.Vertices(I(kk)+1:I(kk+1)-1,:),FloeWhole.SubFloes(ii).h);
        end
    else
        inertia(ii) = PolygonMoments(FloeWhole.SubFloes(ii).poly.Vertices,FloeWhole.SubFloes(ii).h);
    end
    [Xi,Yi] = centroid(FloeWhole.SubFloes(ii).poly);
    centers(ii,:) = [Xi,Yi];
end
FloeWhole.Xm = sum(rho_ice*areaS.*cat(1,FloeWhole.SubFloes.h).*centers(:,1))./FloeWhole.mass;
FloeWhole.Ym = sum(rho_ice*areaS.*cat(1,FloeWhole.SubFloes.h).*centers(:,2))./FloeWhole.mass;
FloeWhole.inertia_moment = sum(inertia+cat(1,FloeWhole.SubFloes.h).*sqrt((centers(:,1)-FloeWhole.Xm).^2+(centers(:,2)-FloeWhole.Ym).^2));

%Calculate values for floe losing ice
[poly2new] = subtract(FloeLess.poly,FloeWhole.poly);
polyout = sortregions(poly2new,'area','descend');
R = regions(polyout);
R = R(area(R)>1e4);
FloeNEW = [];
for kk = 1:length(R)
    FloeLess = Floe2;
    poly2new = R(kk);
    FloeLess.poly = poly2new;
    FloeLess.area = area(poly2new);
    [Xi,Yi] = centroid(poly2new);
    FloeLess.Xi = Xi;
    FloeLess.Yi = Yi;
    FloeLess.mass = FloeLess.mass+V*rho_ice;
    alive = ones(1,length(FloeLess.SubFloes));
    for ii = 1:length(overlap2)
        polynew = intersect(FloeLess.poly,FloeLess.SubFloes(overlap2(ii)).poly);
        [poly2new] = subtract(polynew,FloeWhole.poly);
        if area(poly2new) < 10
            alive(overlap2(ii)) = 0;
        else
            polyout = sortregions(poly2new,'area','descend');
            R2 = regions(polyout);
            poly2new = R2(1);
            poly2new = rmholes(poly2new);
            FloeLess.SubFloes(overlap2(ii)).poly = poly2new;
            if length(R2) > 1
                poly2new = R2(2);
                poly2new = rmholes(poly2new);
                if area(poly2new) > 1000
                    FloeLess.SubFloes(length(FloeLess.SubFloes)+1).poly = poly2new;
                    FloeLess.SubFloes(length(FloeLess.SubFloes)).h = FloeLess.SubFloes(overlap2(ii)).h;
                end
            end
        end
    end
    FloeLess.SubFloes(alive == 0) = [];
    NL = length(FloeLess.SubFloes);
    if NL == 0
        FloeLess.alive = 0;
    end
    areaS = zeros(NL,1);
    inertia = zeros(NL,1);
    centers = zeros(NL,2);
    for ii = 1:NL
        areaS(ii) = area(FloeLess.SubFloes(ii).poly);
        inertia(ii) = PolygonMoments(FloeLess.SubFloes(ii).poly.Vertices,FloeLess.SubFloes(ii).h);
        if FloeLess.SubFloes(ii).poly.NumHoles > 0
            breaks = isnan(FloeLess.SubFloes(ii).poly.Vertices(:,1));
            I = find(breaks == 1);
            I = [0 I' length(breaks)+1];
            inertia(ii) = 0;
            for kk = length(I) -1
                inertia(ii) = inertia(ii) + PolygonMoments(FloeLess.SubFloes(ii).poly.Vertices(I(kk)+1:I(kk+1)-1,:),FloeLess.SubFloes(ii).h);
            end
        else
            inertia(ii) = PolygonMoments(FloeLess.SubFloes(ii).poly.Vertices,FloeLess.SubFloes(ii).h);
        end
        [Xi,Yi] = centroid(FloeLess.SubFloes(ii).poly);
        centers(ii,:) = [Xi,Yi];
    end
    FloeLess.mass = sum(rho_ice*areaS.*cat(1,FloeLess.SubFloes.h));
    FloeLess.Xm = sum(rho_ice*areaS.*cat(1,FloeLess.SubFloes.h).*centers(:,1))./FloeLess.mass;
    FloeLess.Ym = sum(rho_ice*areaS.*cat(1,FloeLess.SubFloes.h).*centers(:,2))./FloeLess.mass;
    FloeLess.inertia_moment = sum(inertia+cat(1,FloeLess.SubFloes.h).*sqrt((centers(:,1)-FloeLess.Xm).^2+(centers(:,2)-FloeLess.Ym).^2));
    FloeLess.rmax = sqrt(max(sum((FloeLess.poly.Vertices' - [FloeLess.Xi; FloeLess.Yi]).^2,1)));

    FloeNEW = [FloeNEW FloeLess];
end

end

