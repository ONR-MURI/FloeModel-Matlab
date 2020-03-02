function [Floe1,Floe2,dissolvedNEW]= ridging(dissolvedNEW,Floe1,Floe2,Nx,Ny,c2_boundary)
%% 
polyout = intersect(Floe1.poly,Floe2.poly);
areaPoly = area(polyout);
rho_ice=920;
rho_l = 997;
E = max([Floe1.E, Floe2.E]);
sigma_m = max([Floe1.sigma_m, Floe2.sigma_m]);
nu = 0.29;
g = 9.81;
hc = 14.2*(1-nu^2)/(rho_l*g)*sigma_m^2/E;
disolved = 0;

%check to make sure one floe is not inside the other
if areaPoly/area(Floe1.poly)>0.9
    dissolvedNEW = dissolvedNEW+calc_vol_dissolved(Floe1,Nx,Ny,c2_boundary);
    disolved = 1;
    Floe1.alive = 0;
elseif areaPoly/area(Floe2.poly)>0.9
    dissolvedNEW = dissolvedNEW+calc_vol_dissolved(Floe2,Nx,Ny,c2_boundary);
    disolved = 1;
    Floe2.alive = 0;
end

if disolved == 0 && areaPoly > 500
    if Floe1.h>= hc && Floe2.h >= hc
        p=1/(1+Floe1.h/Floe2.h);
        if rand(1)>= p
            overlapV = areaPoly*Floe2.h;
            Floe1.h = Floe1.h+overlapV/Floe1.area;
            Floe1.mass = area(Floe1.poly)*Floe1.h*rho_ice;
            Floe1.inertia_moment = PolygonMoments(Floe1.poly.Vertices,Floe1.h);
            [poly2new] = subtract(Floe2.poly,Floe1.poly);
%             if poly2new.NumRegions > 1
                polyout = sortregions(poly2new,'area','descend');
                R = regions(polyout);
                poly2new = R(1);
%             end
            poly2new = translate(poly2new, [(-1)^randi([0 1])*rand (-1)^randi([0 1])*rand]);
            [Floe2.Xi,Floe2.Yi] = centroid(poly2new);
            Floe2.Xi = Floe2.Xi;
            Floe2.Yi = Floe2.Yi;
            Floe2.area = area(poly2new);
            Floe2.mass = area(poly2new)*Floe2.h*rho_ice;
            Floe2.poly = poly2new;
            Floe2.rmax = sqrt(max(sum((poly2new.Vertices' - [Floe2.Xi; Floe2.Yi]).^2,1)));
            Floe2.inertia_moment=PolygonMoments(Floe2.poly.Vertices,Floe2.h);
        else
            overlapV = areaPoly*Floe1.h;
            Floe2.h = Floe2.h+overlapV/Floe2.area;
            Floe2.mass = area(Floe2.poly)*Floe2.h*rho_ice;
            Floe2.inertia_moment=PolygonMoments(Floe2.poly.Vertices,Floe2.h);
            [poly1new] = subtract(Floe1.poly,Floe2.poly);
%             if poly1new.NumRegions > 1
                polyout = sortregions(poly1new,'area','descend');
                R = regions(polyout);
                poly1new = R(1);
%             end
            poly1new = translate(poly1new, [(-1)^randi([0 1])*rand (-1)^randi([0 1])*rand]);
            [Floe1.Xi,Floe1.Yi] = centroid(poly1new);
            Floe1.Xi = Floe1.Xi;
            Floe1.Yi = Floe1.Yi;
            Floe1.area = area(poly1new);
            Floe1.mass = area(poly1new)*Floe1.h*rho_ice;
            Floe1.poly = poly1new;
            Floe1.rmax = sqrt(max(sum((poly1new.Vertices' - [Floe1.Xi; Floe1.Yi]).^2,1)));
            Floe1.inertia_moment = PolygonMoments(Floe1.poly.Vertices,Floe1.h);
        end
    elseif Floe1.h>= hc && Floe2.h< hc
        overlapV = areaPoly*Floe2.h;
        Floe1.h = Floe1.h+overlapV/Floe1.area;
        Floe1.mass = area(Floe1.poly)*Floe1.h*rho_ice;
        Floe1.inertia_moment = PolygonMoments(Floe1.poly.Vertices,Floe1.h);
        [poly2new] = subtract(Floe2.poly,Floe1.poly);
%         if poly2new.NumRegions > 1
            polyout = sortregions(poly2new,'area','descend');
            R = regions(polyout);
            poly2new = R(1);
%         end
        poly2new = translate(poly2new, [(-1)^randi([0 1])*rand (-1)^randi([0 1])*rand]);
        [Floe2.Xi,Floe2.Yi] = centroid(poly2new);
        Floe2.Xi = Floe2.Xi;
        Floe2.Yi = Floe2.Yi;
        Floe2.area = area(poly2new);
        Floe2.mass = area(poly2new)*Floe2.h*rho_ice;
        Floe2.poly = poly2new;
        Floe2.rmax = sqrt(max(sum((poly2new.Vertices' - [Floe2.Xi; Floe2.Yi]).^2,1)));
        Floe2.inertia_moment=PolygonMoments(Floe2.poly.Vertices,Floe2.h);
    elseif Floe1.h < hc && Floe2.h >= hc
        overlapV = areaPoly*Floe1.h;
        Floe2.h = Floe2.h+overlapV/Floe2.area;
        Floe2.mass = area(Floe2.poly)*Floe2.h*rho_ice;
        Floe2.inertia_moment=PolygonMoments(Floe2.poly.Vertices,Floe2.h);
        [poly1new] = subtract(Floe1.poly,Floe2.poly);
%         if poly1new.NumRegions > 1
            polyout = sortregions(poly1new,'area','descend');
            R = regions(polyout);
            poly1new = R(1);
%         end
        poly1new = translate(poly1new, [(-1)^randi([0 1])*rand (-1)^randi([0 1])*rand]);
        [Floe1.Xi,Floe1.Yi] = centroid(poly1new);
        Floe1.Xi = Floe1.Xi;
        Floe1.Yi = Floe1.Yi;
        Floe1.area = area(poly1new);
        Floe1.mass = area(poly1new)*Floe1.h*rho_ice;
        Floe1.poly = poly1new;
        Floe1.rmax = sqrt(max(sum((poly1new.Vertices' - [Floe1.Xi; Floe1.Yi]).^2,1)));
        Floe1.inertia_moment = PolygonMoments(Floe1.poly.Vertices,Floe1.h);
    end
end

if length(Floe1.poly.Vertices) > 500
    Floe1 = FloeSimplify2(Floe1, 250);
    polyout = sortregions(Floe1.poly,'area','descend');
    R = regions(polyout);
    poly1new = R(1);
    poly1new = rmholes(poly1new);
    [Floe1.Xi,Floe1.Yi] = centroid(poly1new);
    Floe1.Xi = Floe1.Xi;
    Floe1.Yi = Floe1.Yi;
    Floe1.area = area(poly1new);
    Floe1.mass = area(poly1new)*Floe1.h*rho_ice;
    Floe1.poly = poly1new;
    Floe1.rmax = sqrt(max(sum((poly1new.Vertices' - [Floe1.Xi; Floe1.Yi]).^2,1)));
    Floe1.inertia_moment = PolygonMoments(Floe1.poly.Vertices,Floe1.h);
elseif length(Floe2.poly.Vertices) > 500
    Floe2 = FloeSimplify2(Floe2, 250);
    polyout = sortregions(Floe2.poly,'area','descend');
    R = regions(polyout);
    poly2new = R(1);
    poly2new = rmholes(poly2new);
    poly2new = translate(poly2new, [(-1)^randi([0 1])*rand (-1)^randi([0 1])*rand]);
    [Floe2.Xi,Floe2.Yi] = centroid(poly2new);
    Floe2.Xi = Floe2.Xi;
    Floe2.Yi = Floe2.Yi;
    Floe2.area = area(poly2new);
    Floe2.mass = area(poly2new)*Floe2.h*rho_ice;
    Floe2.poly = poly2new;
    Floe2.rmax = sqrt(max(sum((poly2new.Vertices' - [Floe2.Xi; Floe2.Yi]).^2,1)));
    Floe2.inertia_moment=PolygonMoments(Floe2.poly.Vertices,Floe2.h);
end

end