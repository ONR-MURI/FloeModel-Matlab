function [floenew] = FloeSimplify2(floe, ddx)
%Take polyshape with a lot of vertices and simplify it to have fewer
%vertices
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
floenew = floe;
x = floe.Xi-floe.rmax:ddx:floe.Xi+floe.rmax;
y = floe.Yi-floe.rmax:ddx:floe.Yi+floe.rmax;
[xx,yy]=meshgrid(x,y);
P = [xx(:) yy(:)];

[k,~] = dsearchn(P,floe.poly.Vertices);
polynew = polyshape(P(k,1),P(k,2));
floenew.poly = simplify(polynew);
warning('on',id)

Volume = floe.area*floe.h;
floenew.area = area(floenew.poly);
floenew.h = Volume/floenew.area;
floenew.rmax = sqrt(max(sum((floenew.poly.Vertices' - [floenew.Xi; floenew.Yi]).^2,1)));
floenew.inertia_moment=PolygonMoments(floenew.poly.Vertices,floenew.h);
end

