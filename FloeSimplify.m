function [floenew] = FloeSimplify(floe, ddx,SUBFLOES)
%Take polyshape with a lot of vertices and simplify it to have fewer
%vertices
%% Remap the main polygon to a shape with fewer vertices
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
rho_ice = 920;
floenew = floe;
n=(fix(floe.rmax/ddx)+1); n=ddx*(-n:n);
[Xg, Yg]= meshgrid(n+floe.Xi, n+floe.Yi);
P = [Xg(:) Yg(:)];

[k,~] = dsearchn(P,floe.poly.Vertices);
polynew = polyshape(P(k,1),P(k,2));
polyout = sortregions(polynew,'area','descend');
R = regions(polyout);
poly1new = R(1);
poly1new = rmholes(poly1new);
floenew.poly = simplify(poly1new);
floenew.area = area(floenew.poly);
[Xi,Yi] = centroid(floenew.poly);
floenew.Xi = Xi; 
floenew.Yi = Yi;

%% Now take care of the subfloes

%Use exisiting points to generate new Voronoi shapes
X = floe.vorX; Y = floe.vorY;
boundingbox = floe.vorbox;
[~, b] = polybnd_voronoi([X Y],boundingbox);
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
[floex, floey] = centroid([floe.SubFloes.poly]);
for ii = 1:N
    [Xi,Yi] = centroid(subfloes(ii));
    [I,~] = dsearchn([floex',floey'],[Xi,Yi]);
    centers(ii,:) = [Xi,Yi];
%     poly = intersect(subfloes(ii),[floe.SubFloes.poly]);
%     [~,I] = min(abs(area(poly)/area(subfloes(ii))-1));
    floenew.SubFloes(ii).poly = subfloes(ii);
    floenew.SubFloes(ii).h = floe.SubFloes(I).h;
    areaS(ii) = area(floenew.SubFloes(ii).poly);
    inertia(ii) = PolygonMoments(floenew.SubFloes(ii).poly.Vertices,floenew.SubFloes(ii).h);
end

%% 

floenew.Xm = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h).*centers(:,1))./floenew.mass;
floenew.Ym = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h).*centers(:,2))./floenew.mass;
floenew.inertia_moment = sum(inertia+cat(1,floenew.SubFloes.h).*sqrt((centers(:,1)-floenew.Xm).^2+(centers(:,2)-floenew.Ym).^2));
floenew.rmax = sqrt(max(sum((floenew.poly.Vertices' - [floenew.Xi; floenew.Yi]).^2,1)));

% if sum(floenew.h.*areaS)/floenew.area<sum([floe.SubFloes.h].*area([floe.SubFloes.poly]))/floe.area && length(floenew.SubFloes)<= length(floe.SubFloes)
%     xx = 1;
%     xx(1) = [1 2];
% end

warning('on',id)
end

