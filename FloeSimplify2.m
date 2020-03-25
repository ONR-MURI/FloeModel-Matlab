function [floenew] = FloeSimplify2(floe, ddx)
%Take polyshape with a lot of vertices and simplify it to have fewer
%vertices
%% 
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
rho_ice = 920;
floenew = floe;
n=(fix(floe.rmax/ddx)+1); n=ddx*(-n:n);
[Xg, Yg]= meshgrid(n+floe.Xi, n+floe.Yi);
% x = floe.Xi-floe.rmax:ddx:floe.Xi+floe.rmax;
% y = floe.Yi-floe.rmax:ddx:floe.Yi+floe.rmax;
% [xx,yy]=meshgrid(x,y);
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
polysmall = scale(floenew.poly, 0.98);
[xs,ys] = centroid(polysmall);
polysmall = translate(polysmall, [floenew.Xi-xs, floenew.Yi-ys]);
polysmall = intersect(polysmall,floenew.poly);
%% 

N = length(floenew.SubFloes);
alive = ones(1,N);
for ii = 1:N
    Volume = area(floe.SubFloes(ii).poly)*floe.SubFloes(ii).h;
    [in, ~] = inpolygon(floe.SubFloes(ii).poly.Vertices(:,1),floe.SubFloes(ii).poly.Vertices(:,2),polysmall.Vertices(:,1),polysmall.Vertices(:,2));
    on = logical(abs(in-1));
    [k,~] = dsearchn(P,floe.SubFloes(ii).poly.Vertices);  
    vertices = floe.SubFloes(ii).poly.Vertices;
    vertices(on,1) = P(k(on),1);
    vertices(on,2) = P(k(on),2);
    polynew = polyshape(vertices);
    polynew = union(polynew, floe.SubFloes(ii).poly);
    polynew = intersect(polynew,floenew.poly);
    if ii == 1
        polyU = polynew;
    else
        polynew = subtract(polynew, polyU);
        polyU = union(polynew,polyU);
    end
    if area(polynew) > 100
        polyout = sortregions(polynew,'area','descend');
        R = regions(polyout);
        poly1 = R(1);
        poly1 = rmholes(poly1);
        floenew.SubFloes(ii).poly =poly1;
        floenew.SubFloes(ii).h = Volume/area(floenew.SubFloes(ii).poly);
        if floenew.SubFloes(ii).h > 30
            floenew.SubFloes(ii).h = 30;
        elseif floenew.SubFloes(ii).h < floe.SubFloes(ii).h
            floenew.SubFloes(ii).h = floe.SubFloes(ii).h;
        end
    else
        alive(ii) = 0;
    end
end
floenew.SubFloes(alive ==0) = [];
warning('on',id)

N1 = length(floenew.SubFloes);
areaS = zeros(N1,1);
inertia = zeros(N1,1);
centers = zeros(N1,2);
for ii = 1:N1
    areaS(ii) = area(floenew.SubFloes(ii).poly);
    inertia(ii) = PolygonMoments(floenew.SubFloes(ii).poly.Vertices,floenew.SubFloes(ii).h);
    [Xi,Yi] = centroid(floenew.SubFloes(ii).poly);
    centers(ii,:) = [Xi,Yi];
end
floenew.Xm = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h).*centers(:,1))./floenew.mass;
floenew.Ym = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h).*centers(:,2))./floenew.mass;
floenew.inertia_moment = sum(inertia+cat(1,floenew.SubFloes.h).*sqrt((centers(:,1)-floenew.Xm).^2+(centers(:,2)-floenew.Ym).^2));
floenew.rmax = sqrt(max(sum((floenew.poly.Vertices' - [floenew.Xi; floenew.Yi]).^2,1)));
end

