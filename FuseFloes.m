function [floenew] = FuseFloes(floe1,floe2,SUBFLOES)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
rho_ice = 920;
height.mean = 1.5;
height.delta = 0;

polynew = union(floe1.poly, floe2.poly);
polyout = sortregions(polynew,'area','descend');
R = regions(polyout);
poly1new = R(1);
poly1new = rmholes(poly1new);
poly = simplify(poly1new);
floenew = initialize_floe_values(poly, height, SUBFLOES);
floenew.area = area(floenew.poly);
[Xi,Yi] = centroid(floenew.poly);
floenew.Xi = Xi; 
floenew.Yi = Yi;

%% Now take care of the subfloes

%Use exisiting points to generate new Voronoi shapes
X = [floe1.vorX; floe2.vorX]; Y = [floe1.vorY; floe2.vorY];
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
    floenew.SubFloes(ii).h = h(I);
    areaS(ii) = area(floenew.SubFloes(ii).poly);
    inertia(ii) = PolygonMoments(floenew.SubFloes(ii).poly.Vertices,floenew.SubFloes(ii).h);
end

%% 
floenew.mass = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h));
floenew.Xm = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h).*centers(:,1))./floenew.mass;
floenew.Ym = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h).*centers(:,2))./floenew.mass;
floenew.inertia_moment = sum(inertia+cat(1,floenew.SubFloes.h).*sqrt((centers(:,1)-floenew.Xm).^2+(centers(:,2)-floenew.Ym).^2));
floenew.rmax = sqrt(max(sum((floenew.poly.Vertices' - [floenew.Xi; floenew.Yi]).^2,1)));

floenew.Ui = (floe1.Ui*floe1.mass + floe2.Ui*floe2.mass)/(floe1.mass+floe2.mass);
floenew.Vi = (floe1.Vi*floe1.mass + floe2.Vi*floe2.mass)/(floe1.mass+floe2.mass);
floenew.dUi_p = (floe1.dUi_p*floe1.mass + floe2.dUi_p*floe2.mass)/(floe1.mass+floe2.mass);
floenew.dVi_p = (floe1.dVi_p*floe1.mass + floe2.dVi_p*floe2.mass)/(floe1.mass+floe2.mass);
floenew.ksi_ice = (floe1.ksi_ice*floe1.mass + floe2.ksi_ice*floe2.mass)/(floe1.mass+floe2.mass);
floenew.potentialInteractions = floe1.potentialInteractions;

end

