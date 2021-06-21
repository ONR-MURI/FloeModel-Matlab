function [Floes] = frac_corner(floe,grind,poly)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
% id3 ='MATLAB:polyshape:boundary3Points';
% warning('off',id3)


rho_ice = 920;
Floes = [];
fracturedFloes = [];
angles = polyangles(poly.Vertices(:,1),poly.Vertices(:,2));
Anorm = 180-360/length(angles);
angs = unique(angles(grind));
if length(angs) < length(angles(grind)) && grind(1)==1 && grind(end) == 1
    grind(end) = 0;
end
alph = min(angles);
for jj = 1:length(grind)
    if grind(jj)
        ii = jj;
        if ii>1 && ii<length(poly.Vertices)
            X = poly.Vertices(ii-1:ii+1,1); Y = poly.Vertices(ii-1:ii+1,2);
        elseif ii == 1
            X = [poly.Vertices(end,1); poly.Vertices(1:2,1)];
            Y = [poly.Vertices(end,2); poly.Vertices(1:2,2)];
        else
            X = [poly.Vertices(end-1:end,1); poly.Vertices(1,1)];
            Y = [poly.Vertices(end-1:end,2); poly.Vertices(1,2)];
        end
        d1 = sqrt((X(1)-X(2))^2+(Y(1)-Y(2))^2); d2 = sqrt((X(3)-X(2))^2+(Y(3)-Y(2))^2);
        d = min(d1,d2);
%         xx = 1; xx(1) =[1 2];
        if abs(X(1)-X(2)) < 120 && abs(Y(1)-Y(2)) < 120
            x1 = X(1); y1 = Y(1);
        elseif abs(d*alph/Anorm/5) < 120 && abs(d*alph/Anorm/5)<120
            x1 = X(2)+120/d1*(X(1)-X(2)); y1 = Y(2)+120/d1*(Y(1)-Y(2));
        else
            x1 = X(2)+d/d1*alph/Anorm*(X(1)-X(2))/5; y1 = Y(2)+d/d1*alph/Anorm*(Y(1)-Y(2))/5;
        end
        if abs(X(3)-X(2)) < 120 && abs(Y(3)-Y(2)) < 120
            x2 = X(3); y2 = Y(3);
        elseif abs(d*alph/Anorm/5) < 120 && abs(d*alph/Anorm/5)<120
            x2 = X(2)+120/d2*(X(3)-X(2)); y2 = Y(2)+120/d2*(Y(3)-Y(2));
        else
            x2 = X(2)+d/d2*alph/Anorm*(X(3)-X(2))/5; y2 = Y(2)+d/d2*alph/Anorm*(Y(3)-Y(2))/5;
        end
        if abs(x1-x2)<0.1
            x1 = x1+1;
        end
        if ii>1 && ii<length(poly.Vertices)
            Pc1 = [poly.Vertices(1:ii-1,:); x1 y1; x2 y2; poly.Vertices(ii+1:end,:)];
            Pc2 = [x1 y1; x2 y2; poly.Vertices(ii,:)];
            poly1 = simplify(polyshape(Pc1));
            poly2 = simplify(polyshape(Pc2));
            if poly1.NumRegions > 1
                Pc1 = [poly.Vertices(1:ii-1,:); x2 y2; x1 y1; poly.Vertices(ii+1:end,:)];
                poly1 = simplify(polyshape(Pc1));
                %         if poly1.NumRegions > 1
                %             xx = 1; xx(1) = [ 1 2];
                %         end
            end
        elseif ii == 1
            Pc1 = [ x1 y1; x2 y2; poly.Vertices(ii+1:end,:)];
            Pc2 = [x1 y1; x2 y2; poly.Vertices(ii,:)];
            poly1 = simplify(polyshape(Pc1));
            poly2 = simplify(polyshape(Pc2));
            if poly1.NumRegions > 1
                Pc1 = [x2 y2; x1 y1; poly.Vertices(ii+1:end,:)];
                poly1 = simplify(polyshape(Pc1));
                %         if poly1.NumRegions > 1
                %             xx = 1; xx(1) = [ 1 2];
                %         end
            end
        else
            Pc1 = [poly.Vertices(1:ii-1,:); x1 y1; x2 y2];
            Pc2 = [x1 y1; x2 y2; poly.Vertices(ii,:)];
            poly1 = simplify(polyshape(Pc1));
            poly2 = simplify(polyshape(Pc2));
            if poly1.NumRegions > 1
                Pc1 = [poly.Vertices(1:ii-1,:); x2 y2; x1 y1];
                poly1 = simplify(polyshape(Pc1));
                %         if poly1.NumRegions > 1
                %             xx = 1; xx(1) = [ 1 2];
                %         end
            end
        end
        % Pc1 = cutpolygon([poly.Vertices;poly.Vertices(1,:)],[x1, y1; x2, y2],1);
        % Pc2 = cutpolygon([poly.Vertices;poly.Vertices(1,:)],[x1, y1; x2, y2],2);
        % poly1 = simplify(polyshape(Pc1));
        % poly2 = simplify(polyshape(Pc2));
        % xx = 1; xx(1) = [1 2];
        
        R1 = [regions(poly1); regions(poly2)];
        a = R1(area(R1)>10);
        Atot = sum(area(a));
        Mtot = floe.mass*Atot/area(polyshape(floe.c_alpha'));
        a = poly2;
        
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
            
            FloeNEW.Xi = floe.Xi+Xi; FloeNEW.Yi = floe.Yi+Yi; FloeNEW.alive = 1;
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
            
            fracturedFloes = [fracturedFloes FloeNEW];
            clear FloeNEW
        end
%         xx = 1; xx(1) =[1 2];
    end
end

if ~isempty(fracturedFloes)
    polyu = union([fracturedFloes.poly]); polyu = translate(polyu,-[floe.Xi, floe.Yi]);
    polynew = subtract(poly,polyu);
    if polynew.NumRegions > 1
        polyout = sortregions(polynew,'area','descend');
        R = regions(polyout);
        polynew = R(1);
    end
    FloeNEW.poly = rmholes(polynew);
    [Xi,Yi] = centroid(FloeNEW.poly);
    FloeNEW.area = area(FloeNEW.poly);
    FloeNEW.mass = floe.mass*area(polynew)/area(poly);
    FloeNEW.h = FloeNEW.mass/(rho_ice*FloeNEW.area);
    FloeNEW.c_alpha = [(FloeNEW.poly.Vertices-[Xi Yi])' [FloeNEW.poly.Vertices(1,1)-Xi; FloeNEW.poly.Vertices(1,2)-Yi]];
    FloeNEW.c0 = FloeNEW.c_alpha;
    FloeNEW.inertia_moment = PolygonMoments(FloeNEW.c0',FloeNEW.h);
    FloeNEW.angles = polyangles(FloeNEW.poly.Vertices(:,1),FloeNEW.poly.Vertices(:,2));
    FloeNEW.rmax = sqrt(max(sum((FloeNEW.poly.Vertices' - [Xi;Yi]).^2,1)));
    FloeNEW.poly = translate(FloeNEW.poly,[floe.Xi,floe.Yi]);
    FloeNEW.strain=floe.strain;
    FloeNEW.Stress = floe.Stress*area(polynew)/area(poly);
    FloeNEW.MaxShear = floe.MaxShear*area(polynew)/area(poly);
    FloeNEW.Fx = floe.Fx; FloeNEW.Fy = floe.Fy;
    FloeNEW.FxOA = 0; FloeNEW.FyOA = 0; FloeNEW.torqueOA = 0;
    
    FloeNEW.X = FloeNEW.rmax*(2*rand(1000,1) - 1);
    FloeNEW.Y = FloeNEW.rmax*(2*rand(1000,1) - 1);
    FloeNEW.A = inpolygon(FloeNEW.X,FloeNEW.Y,FloeNEW.c_alpha(1,:),FloeNEW.c_alpha(2,:));
    
    FloeNEW.Xi = floe.Xi+Xi; FloeNEW.Yi = floe.Yi+Yi; FloeNEW.alive = 1;
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
    FloeNEW.collision_torque = 0;
    FloeNEW.OverlapArea = 0;
end

Floes = [FloeNEW fracturedFloes];

floenew = [];
for ii = 1:length(Floes)
    floe2 = FloeSimplify(Floes(ii),false);
    for jj = 1:length(floe2)
        if jj == 1
            Floes(ii) = floe2(jj);
        else
            floenew = [floenew floe2(jj)];
        end
    end
end
Floes = [Floes floenew];

% warning('on',id)
% warning('on',id3)
% if floe.mass/sum(cat(1,Floes.mass))-1 > 1e-3
%     xx = 1;
%     xx(1) = [1 2];
% end
% 
% if sum(cat(1,Floes.mass))/floe.mass-1 > 1e-3
%     xx = 1;
%     xx(1) = [1 2];
% end
% 
% h = cat(1,Floes.h);
% if max(h)/floe.h-1 > 1e-2
%     xx = 1;
%     xx(1) = [1 2];
% end
end

