function [Floe1,Floe2]= ridging(Floe1,Floe2,c2_boundary_poly,PERIODIC,min_floe_size)
%% This function takes in two floes and based upon the thickness of the two floes will perform a ridging operation
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)
%Create polyshapes for the boundary in a nonperiodic case as well as polyshapes for the union of all subfloes 
floe1 = Floe1;
floe2 = Floe2;

poly1 = simplify(polyshape(Floe1.c_alpha'+[Floe1.Xi Floe1.Yi]));
poly2 = simplify(polyshape(Floe2.c_alpha'+[Floe2.Xi Floe2.Yi]));
Floe1.poly = poly1;
Floe2.poly = poly2;


%Find area of overlap
aPoly = area(intersect(poly1,poly2));

%Find critical thickness
rho_ice=920;
rho_l = 997;
E = 1e9;
sigma_m = 400000;
nu = 0.29;
g = 9.81;
hc = 0.5;%14.2*(1-nu^2)/(rho_l*g)*sigma_m^2/E;
disolved = 0;
boundary = 0;

%check to make sure one floe is not inside the other and add mass to
%dissolved if it is
if area(poly2)/area(c2_boundary_poly)>0.99
    boundary = 1;
elseif aPoly/area(Floe1.poly)>0.75 || Floe1.area<min_floe_size
    disolved = 1;
    Floe1.alive = 0;
elseif aPoly/area(Floe2.poly)>0.75 || Floe2.area < min_floe_size
    disolved = 1;
    Floe2.alive = 0;
elseif Floe1.alive+Floe2.alive < 2
    disolved = 0;
end
%% If there is enough overlap then allow ridging to happen
if ~isfield(Floe2,'poly')
    xx = 1;
    xx(1) = [1 2];
end
if disolved == 0 && aPoly > 500 && ~boundary
    
    %Determine overlap in the different subfloes
    V1 = aPoly*Floe1.h;
    V2 = aPoly*Floe2.h;
    
    %Use the thicknesses to determine how mass will be transfered
    if Floe1.h>= hc && Floe2.h >= hc
        p=1/(1+Floe1.h/Floe2.h);
        if rand(1)>= p
            [Floe1, Floe2] = ridge_values_update(Floe1,Floe2, V2);
        else
            [Floe2, Floe1] = ridge_values_update(Floe2,Floe1, V1);
        end
    elseif Floe1.h>= hc && Floe2.h< hc
        [Floe1, Floe2] = ridge_values_update(Floe1,Floe2, V2);
    elseif Floe1.h < hc && Floe2.h >= hc
        [Floe2, Floe1] = ridge_values_update(Floe2,Floe1, V1);
    end
end
if isempty(Floe2)
    
elseif ~isfield(Floe2,'poly')
    xx = 1;
    xx(1) = [1 2];
end

%Perform ridging with boundary
if ~PERIODIC && boundary
    Lx= max(c2_boundary_poly.Vertices(:,1));
    Ly= max(c2_boundary_poly.Vertices(:,2));%c2 must be symmetric around x=0 for channel boundary conditions.
    x=[-1 -1 1 1 -1]*Lx*2;
    y=[-1 1 1 -1 -1]*Ly*2;
    polybound = polyshape(x,y);
    c2_poly = subtract(polybound,c2_boundary_poly);
    Abound = area(intersect(Floe1.poly,c2_poly));
    poly1new = subtract(Floe1.poly,c2_poly);
    V1 = Abound*Floe1.h; 

    if area(poly1new) < 1000
        Floe1.alive = 0;
    elseif Abound>0 && area(poly1new) >= 1000
        %Calculate new floe properties after the mass transfer and shape is
        %updated
        polyout = sortregions(poly1new,'area','descend');
        R = regions(polyout);
        poly1new = R(1);
        Floe1.area = area(poly1new);
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

% if max(cat(1,Floe1.h))/floe1.h-1 > 0.15
%     xx = 1;
%     xx(1) = [1 2];
% elseif max(cat(1,Floe2.h))/floe2.h-1 > 0.15
%     xx = 1;
%     xx(1) = [1 2];
% end

if isfield(Floe1,'poly')
    Floe1=rmfield(Floe1,{'poly'});
end
if isfield(Floe2,'poly')
    Floe2=rmfield(Floe2,{'poly'});
else
    xx = 1;
    xx(1) = [1 2];
end

warning('on',id)
warning('on',id3)
end