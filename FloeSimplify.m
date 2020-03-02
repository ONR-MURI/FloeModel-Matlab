function [floenew] = FloeSimplify(floe, c2_boundary, ddx)
%Take polyshape with a lot of vertices and simplify it to have fewer
%vertices
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
floenew = floe;
x = min(c2_boundary(1,:)):ddx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):ddx:max(c2_boundary(2,:));
BW = poly2mask((floe.poly.Vertices(:,1)+max(x))/ddx,(floe.poly.Vertices(:,2)+max(y))/ddx,length(y), length(x));
stats = regionprops(BW);
floenew.Xi = stats.Centroid(1)*ddx-max(x)-ddx;
floenew.Yi = stats.Centroid(2)*ddx-max(y)-ddx;
Anew = logical(BW);
xx=1:size(Anew,2); yy=1:size(Anew,1);
[xx,yy]=meshgrid(ddx*xx,ddx*yy); %% floePixelSize m is the pixel size from Earthdata NASA webpage
xc=mean(xx(Anew)); yc=mean(yy(Anew));
r=sqrt((xx-xc).^2+(yy-yc).^2); r_max=max(r(Anew))+2*ddx; % pad with extra row/column of zeros around the floe

Ngr=2*fix(r_max/ddx); % number of grid boxes over floe

[XX,YY]=meshgrid((-1:2/Ngr:1)*r_max, (-1:2/Ngr:1)*r_max); % floe grid
Xg=(-1:2/Ngr:1)*r_max; Yg=Xg;

A=interp2(xx-xc,yy-yc,double(Anew),XX,YY); A(isnan(A))=0;
A = logical(A);
c0=contourc(Xg,Yg,double(A),[0.5 0.5]);
polynew = polyshape(c0(1,:)+floenew.Xi,c0(2,:)+floenew.Yi);

floenew.poly = simplify(polynew);
warning('on',id)

Volume = floe.area*floe.h;
floenew.area = area(floenew.poly);
floenew.h = Volume/floenew.area;
floenew.rmax = sqrt(max(sum((floenew.poly.Vertices' - [floenew.Xi; floenew.Yi]).^2,1)));
floenew.inertia_moment=PolygonMoments(floenew.poly.Vertices,floenew.h);
end

