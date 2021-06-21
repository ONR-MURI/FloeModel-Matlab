function [Floe1,Floe2] = rafting(Floe1,Floe2,boundary)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
poly1 = simplify(polyshape(Floe1.c_alpha'+[Floe1.Xi Floe1.Yi]));
poly2 = simplify(polyshape(Floe2.c_alpha'+[Floe2.Xi Floe2.Yi]));
Floe1.poly = poly1; Floe2.poly = poly2;
aPoly = area(intersect(poly1,poly2));

%Determine overlap in the different subfloes
if isempty(aPoly)
    V1=0; V2 = 0;
else
    V1 = aPoly*Floe1.h;
    V2 = aPoly*Floe2.h;
end

%Use the thicknesses to determine how mass will be transfered
if V1>0 && ~boundary
    p=1/(1+Floe1.h/Floe2.h);
    if rand(1)>= p
        [Floe1, Floe2] = ridge_values_update(Floe1,Floe2, V2);
    else
        [Floe2, Floe1] = ridge_values_update(Floe2,Floe1, V1);
    end
end

if V1 > 0 && boundary
    polynew = subtract(poly1,poly2);
    Abound = area(intersect(poly1,poly2));
    V1 = Abound*Floe1.h; 

    if area(polynew) < 1000
        Floe1.alive = 0;
    elseif Abound>0 && area(polynew) >= 1000
        %Calculate new floe properties after the mass transfer and shape is
        %updated
        polyout = sortregions(polynew,'area','descend');
        R = regions(polyout);
        poly1new = rmholes(R(1));
        Floe1.area = area(polynew);
        Floe1.h = Floe1.h+V1/Floe1.area;
        Floe1.poly = poly1new;
        [Xi,Yi] = centroid(poly1new);
        Floe1.Xi = Xi;
        Floe1.Yi = Yi;
        Floe1.dalpha_i_p = 0;
        Floe1.alpha_i = 0;
        Floe1.c_alpha = [(poly1new.Vertices-[Floe1.Xi Floe1.Yi])' [poly1new.Vertices(1,1)-Floe1.Xi; poly1new.Vertices(1,2)-Floe1.Yi]];
        Floe1.c0 = Floe1.c_alpha;
        Floe1.angles = polyangles(poly1new.Vertices(:,1),poly1new.Vertices(:,2));
        if poly1new.NumRegions > 1
            xx = 1;
            xx(1) = [1 2];
        end
        
        Floe1.inertia_moment = PolygonMoments(Floe1.c_alpha',Floe1.h);
        Floe1.rmax = max(sqrt(Floe1.c_alpha(1,:).^2+Floe1.c_alpha(2,:).^2));
        Floe1.X = Floe1.rmax*(2*rand(1000,1) - 1);
        Floe1.Y = Floe1.rmax*(2*rand(1000,1) - 1);
        Floe1.A = inpolygon(Floe1.X,Floe1.Y,Floe1.c_alpha(1,:),Floe1.c_alpha(2,:));
    end
end

Floe1=rmfield(Floe1,{'poly'});
Floe2=rmfield(Floe2,{'poly'});

end
