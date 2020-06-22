function [floenew] = FuseFloes(floe1,floe2,SUBFLOES)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
SHIFT = false;
rho_ice = 920;
height.mean = 1.5;
height.delta = 0;

polynew = union(floe1.poly, floe2.poly);
polyout = sortregions(polynew,'area','descend');
R = regions(polyout);
poly1new = R(1);
poly1new = rmholes(poly1new);
poly = simplify(poly1new);
floenew = initialize_floe_values(poly, height,SHIFT,  SUBFLOES);
floenew.area = area(floenew.poly);
[Xi,Yi] = centroid(floenew.poly);
floenew.Xi = Xi; 
floenew.Yi = Yi;
floenew.mass = floe1.mass + floe2.mass;

%% Now take care of the subfloes

%Use exisiting points to generate new Voronoi shapes
[~,d] = dsearchn([floe1.vorX,floe1.vorY],[floe2.vorX,floe2.vorY]);
X = [floe1.vorX; floe2.vorX(d>0)]; Y = [floe1.vorY; floe2.vorY(d>0)];
Xmin = min([floe1.vorbox(:,1);floe2.vorbox(:,1)]); Xmax = max([floe1.vorbox(:,1);floe2.vorbox(:,1)]);
Ymin = min([floe1.vorbox(:,2);floe2.vorbox(:,2)]); Ymax = max([floe1.vorbox(:,2);floe2.vorbox(:,2)]);
floenew.vorbox = [Xmin, Ymin; Xmin Ymax; Xmax, Ymax; Xmax, Ymin];
[~, b] = polybnd_voronoi([X Y],floenew.vorbox);
k = 1;

if SUBFLOES
    for i=1:length(b)
        if ~isnan(b{i}(1))
            a=regions(intersect(floenew.poly,polyshape(b{i})));
            if ~isempty(area(a))
                for j=1:length(a)
                    subfloes(k) = a(j);
                    k = k+1;
                end
            end
        end
    end
else
    subfloes = floenew.poly;
end

%populate new floes with heights of exisiting floes
N = length(subfloes);
floenew.SubFloes = [];
areaS = zeros(N,1);
inertia = zeros(N,1);
centers = zeros(N,2);
[floex1, floey1] = centroid([floe1.SubFloes.poly]);
[floex2, floey2] = centroid([floe2.SubFloes.poly]);
floex = [floex1 floex2]; floey = [floey1 floey2];
h = [[floe1.SubFloes.h] [floe2.SubFloes.h]];
for ii = 1:N
    [Xi,Yi] = centroid(subfloes(ii));
    [I,~] = dsearchn([floex',floey'],[Xi,Yi]);
    centers(ii,:) = [Xi,Yi];
%     poly = intersect(subfloes(ii),[floe.SubFloes.poly]);
%     [~,I] = min(abs(area(poly)/area(subfloes(ii))-1));
    floenew.SubFloes(ii).poly = subfloes(ii);
    if SUBFLOES
        floenew.SubFloes(ii).h = h(I);
    else
        floenew.SubFloes(ii).h = (h(1)*floe1.mass+h(2)*floe2.mass)/floenew.mass;
    end
    areaS(ii) = area(floenew.SubFloes(ii).poly);
    if floenew.SubFloes(ii).poly.NumHoles > 0
        breaks = isnan(floenew.SubFloes(ii).poly.Vertices(:,1));
        I = find(breaks == 1);
        I = [0 I' length(breaks)+1];
        inertia(ii) = 0;
        for kk = length(I) -1
            inertia(ii) = inertia(ii) + PolygonMoments(floenew.SubFloes(ii).poly.Vertices(I(kk)+1:I(kk+1)-1,:),floenew.SubFloes(ii).h);
        end
    else
        inertia(ii) = PolygonMoments(floenew.SubFloes(ii).poly.Vertices,floenew.SubFloes(ii).h);
    end
end

%% 
% floenew.mass = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h));
% test = abs(floenew.mass-floe1.mass-floe2.mass)/(floe1.mass+floe2.mass);
% if test>0.01
%     x = 1;
%     x(1) = [1 2];
% end
floenew.Xm = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h).*centers(:,1))./floenew.mass;
floenew.Ym = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h).*centers(:,2))./floenew.mass;
floenew.inertia_moment = sum(inertia+cat(1,floenew.SubFloes.h).*sqrt((centers(:,1)-floenew.Xm).^2+(centers(:,2)-floenew.Ym).^2));
floenew.rmax = sqrt(max(sum((floenew.poly.Vertices' - [floenew.Xi; floenew.Yi]).^2,1)));

floenew.Ui = (floe1.Ui*floe1.mass + floe2.Ui*floe2.mass)/(floenew.mass);
floenew.Vi = (floe1.Vi*floe1.mass + floe2.Vi*floe2.mass)/(floenew.mass);
floenew.dUi_p = (floe1.dUi_p*floe1.mass + floe2.dUi_p*floe2.mass)/(floenew.mass);
floenew.dVi_p = (floe1.dVi_p*floe1.mass + floe2.dVi_p*floe2.mass)/(floenew.mass);
floenew.ksi_ice = (floe1.ksi_ice*floe1.inertia_moment + floe2.ksi_ice*floe2.inertia_moment)/(floenew.inertia_moment);%use inertia moment instead of mass
floenew.dXi_p = (floe1.dXi_p*floe1.mass + floe2.dXi_p*floe2.mass)/(floenew.mass);
floenew.dYi_p = (floe1.dYi_p*floe1.mass + floe2.dYi_p*floe2.mass)/(floenew.mass);
floenew.dksi_ice_p = (floe1.dksi_ice_p*floe1.inertia_moment + floe2.dksi_ice_p*floe2.inertia_moment)/(floenew.inertia_moment);
floenew.h = (floe1.h*floe1.mass+floe2.h*floe2.mass)/floenew.mass;%floenew.mass/(floenew.area*rho_ice);
floenew = rmfield(floenew,'potentialInteractions');
% floenew.potentialInteractions = floe1.potentialInteractions;


if floenew.poly.NumRegions>1
    xx = 1;
    xx(1) = [1 2];
end

end

