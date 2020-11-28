function [floes] = FloeSimplify(floe)
%Take polyshape with a lot of vertices and simplify it to have fewer
%vertices
%% Remap the main polygon to a shape with fewer vertices
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id2 = 'MATLAB:polyshape:boolOperationFailed';
warning('off',id2)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)
floes = [];
rho_ice = 920;

%Set function to define the tolerance for Douglas-Peuker algorithm
SimpMin = @(A) log10(A)^3.5;

%Create new simplified polyshape
floenew = floe;
%vertnew = DouglasPeucker(floe.poly.Vertices,SimpMin(floe.area));
[vertx,verty] = reducem(floe.c0(1,:)',floe.c0(2,:)',120);
pnew = polyshape(vertx, verty);
Atot = sum(area(pnew));
if isinf(Atot)
  save('floefail.mat','floe','pnew');
%    1
elseif isnan(Atot)
  save('floefail.mat','floe','pnew');
%   2
elseif Atot ==0
  save('floefail.mat','floe','pnew');
%   3
  R = [];
else
  polynew = scale(pnew,sqrt(floe.area/Atot));

  %Align center of old polygon with the enw one
  [x1,y1] = centroid(polynew);
  dx = floe.Xi-x1;
  dy = floe.Yi-y1;
  % if area(intersect(floe.poly,polynew))/floe.area > 0.95
  %     polynew = translate(polynew,[dx, dy]);
  % end

  %Check if simplifcation led to polygon having multiple regions
  polyout = sortregions(polynew,'area','descend');
  R = regions(polyout);
  R = R(area(R)>1e4);
  Atot = sum(area(R));
end

%% Calculate the new properties associated with this floe since it has a new shape
if length(R) == 1
    floes = floenew;
    R = rmholes(R);
    floes.area = area(R);
    [Xi, Yi] = centroid(R);
    floes.Xi = Xi+floenew.Xi; floes.Yi = Yi+floenew.Yi;
    floes.h = floes.mass/(rho_ice*floes.area);
    floes.c0 = [R.Vertices(:,1),R.Vertices(:,2); R.Vertices(1,1),R.Vertices(1,2)]';
    floes.angles = polyangles(R.Vertices(:,1),R.Vertices(:,2));
    A_rot=[cos(floes.alpha_i) -sin(floes.alpha_i); sin(floes.alpha_i) cos(floes.alpha_i)]; %rotation matrix
    floes.c_alpha=A_rot*floes.c0;
    floes.inertia_moment = PolygonMoments(floes.c_alpha',floes.h);
    floes.rmax = max(sqrt(floes.c_alpha(1,:).^2+floes.c_alpha(2,:).^2));
elseif length(R) > 1
    for ii = 1:length(R)
        floenew = floe;
        poly1new = R(ii);
        poly1new = rmholes(poly1new);
        floenew.poly = poly1new;
        [Xi,Yi] = centroid(poly1new);
        floenew.area = area(poly1new);
        floenew.mass = floe.mass*area(R(ii))/Atot;
        floenew.h = floenew.mass/(rho_ice*floenew.area);
        floenew.inertia_moment = PolygonMoments(poly1new.Vertices,floenew.h);
        floenew.c_alpha = [(floenew.poly.Vertices-[Xi Yi])' [floenew.poly.Vertices(1,1)-Xi; floenew.poly.Vertices(1,2)-Yi]];
        floenew.angles = polyangles(poly1new.Vertices(:,1),poly1new.Vertices(:,2));
        floenew.c0 = floenew.c_alpha;
        floenew.rmax = sqrt(max(sum((poly1new.Vertices' - [Xi;Yi]).^2,1)));
        floenew.Xg = floe.Xg;
        floenew.Yg = floe.Yg;
        floenew.X = floe.X;
        floenew. Y = floe.Y;
        
        [in] = inpolygon(floenew.X(:)+Xi, floenew.Y(:)+Yi,floenew.poly.Vertices(:,1),floenew.poly.Vertices(:,2));
        floenew.A=reshape(in,length(floenew.X),length(floenew.X));
        
        floenew.Xi = floe.Xi+Xi; floenew.Yi = floe.Yi+Yi; floenew.alive = 1;
        floenew.alpha_i = 0; floenew.Ui = floe.Ui; floenew.Vi = floe.Vi;
        floenew.dXi_p = floe.dXi_p; floenew.dYi_p = floe.dYi_p;
        floenew.dUi_p = floe.dUi_p; floenew.dVi_p = floe.dVi_p;
        floenew.dalpha_i_p = 0; floenew.ksi_ice = floenew.area/floe.area*floe.ksi_ice;
        floenew.dksi_ice_p = floe.dksi_ice_p;
        floenew.interactions = [];
        floenew.potentialInteractions = [];
        floenew.collision_force = 0;
        floenew.collision_torque = 0;
        floenew.OverlapArea = 0;
                    
        floes = [floes floenew];
    end
else
    floes = floe;
    floes.alive = 0;
end

if floe.mass/sum(cat(1,floes.mass))-1 > 1e-3
    xx = 1;
    xx(1) = [1 2];
end
if sum(cat(1,floes.mass))/floe.mass-1 > 1e-3
    xx = 1;
    xx(1) = [1 2];
end
h = cat(1,floes.h);
if max(h)/floe.h-1 > 0.01
    xx = 1;
    xx(1) = [1 2];
end
% for ii = 1:length(floes)
%     if isempty(floes(ii).SubFloes.inertia)
%         xx=1;
%         xx(1) = [1 2];
%     end
% end
if abs(floe.area/sum(cat(1,floes.area))-1)>0.05
    xx = 1;
    xx(1) =[1 2];
end
for ii = 1:length(floes)
    if abs(floes(ii).area/area(polyshape(floes(ii).c_alpha'))-1)>1e-3
        xx = 1;
        xx(1) =[1 2];
    end
end

warning('on',id)
warning('on',id2)
warning('on',id3)
end

