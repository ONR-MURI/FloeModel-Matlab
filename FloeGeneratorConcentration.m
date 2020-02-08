function [Floe2]= FloeGeneratorConcentration(Floe,c2_boundary,Target,N)
%Generate set of random points
ddx = 250;
id ='MATLAB:polyshape:boundary3Points';
warning('off',id)
%load FloeVoronoi;
N = floor(4*N);%floor(((max(c2_boundary(2,:))-min(c2_boundary(2,:)))*(max(c2_boundary(2,:))-min(c2_boundary(2,:))))/(mean(cat(1,Floe.area))/4));
x = min(c2_boundary(1,:)):ddx:max(c2_boundary(1,:));
y = min(c2_boundary(2,:)):ddx:max(c2_boundary(2,:));
dx = max(x)-min(x);
dy = max(y)-min(y);
X = min(c2_boundary(1,:))-dx/2+2*rand(1,N)*dx;
Y = min(c2_boundary(2,:))-dy/2+2*rand(1,N)*dy;

%Remove any points form inside current floes
for ii = 1:length(Floe)
%     edge = [1:length(Floe(ii).poly.Vertices); 2:length(Floe(ii).poly.Vertices) 1];
%     [stat,~] = inpoly2([X',Y'],Floe(ii).poly.Vertices,edge') ;
%     X(stat) = [];
%     Y(stat) = [];
    [In] =inpolygon(X,Y,Floe(ii).poly.Vertices(:,1),Floe(ii).poly.Vertices(:,2));
    X(In) = [];
    Y(In) = [];
end

%Generate new polygons with voronoi and trim vertices so that they are
%within the boundary
dt = delaunayTriangulation(X',Y');
[Vold,C] = voronoiDiagram(dt);
V = Vold;
[vxold,vyold]=voronoi(X,Y);
[vx, vy] = vertexcrop(x,y,vxold,vyold);
for ii = 1:length(vx)
    [M(ii),jj] = min(abs(V(:,1)-vxold(1,ii)));
    if abs(M(ii))==0
        V(jj,1)=vx(1,ii);
        V(jj,2)=vy(1,ii);
    end
end
poly2 = polyshape(c2_boundary(1,:),c2_boundary(2,:));
xb = [min(x)-dx min(x)-dx max(x)+dx max(x)+dx min(x)-dx];
yb = [min(y)-dy max(y)+dy max(y)+dy min(y)-dy min(y)-dy];
poly3 = polyshape(xb,yb);
polyout = subtract(poly3,poly2);

%calculate current concentration
live = cat(1,Floe.alive);
Floe(logical(abs(live-1)))= [];
Area = sum(cat(1,Floe.area));
cnow = Area/((max(c2_boundary(2,:))-min(c2_boundary(2,:)))*(max(c2_boundary(1,:))-min(c2_boundary(1,:))));

%randomize order for polygons to be added as floes
R = randperm(length(C));
% R = randperm(length(FloeV));
%% Loop until target concentration is reached
Anew = 0;
ii = 1;
ind = 1;
x=cat(1,Floe.Xi);
y=cat(1,Floe.Yi);
rmax=cat(1,Floe.rmax);
alive=cat(1,Floe.alive);
while cnow < Target
    %Take set of vertices from Voronoi and generate polyshape for floe
    if ii>length(R)
        Anew = 0;
    elseif R(ii) > length(C)
%     elseif R(ii) > length(FloeV)
        Anew = 0;
    else
        if C{R(ii)}(1) ==1 %avoid V(1,:) which is inf
            C{R(ii)}(1)=[];
            C{R(ii)} = C{R(ii)}(2:end);
        end
        Vnow1 = V(C{R(ii)},1);
        Vnow1(isnan(Vnow1))=[];
        Vnow2 = V(C{R(ii)},2);
        Vnow2(isnan(Vnow2))=[];
        poly1 = polyshape(Vnow1,Vnow2);
%         poly1 = FloeV(R(ii));
        
        % Generator Potential Interactions for newly created floes
        FloeNEW(ii).potentialInteractions=[];
        [xi, yi] = centroid(poly1);
        r_max = sqrt(max(sum((poly1.Vertices' - [xi;yi]).^2,1)));
        k=1;
        if area(poly1)>0
            for j=1:length(Floe)
                if  alive(j) && sqrt((xi-x(j))^2 + (yi-y(j))^2)<(r_max+rmax(j)) % if floes are potentially overlapping
                    if area(poly1)>0
                        FloeNEW(ii).potentialInteractions(k).floeNum=j;
                        FloeNEW(ii).potentialInteractions(k).c=Floe(j).poly;
%                         [In] =inpolygon(poly1.Vertices(:,1),poly1.Vertices(:,2),Floe(jj).poly.Vertices(:,1),Floe(jj).poly.Vertices(:,2));
                        edge = [1:length(Floe(j).poly.Vertices); 2:length(Floe(j).poly.Vertices) 1];
                        [stat,~] = inpoly2([poly1.Vertices(:,1),poly1.Vertices(:,2)],Floe(j).poly.Vertices,edge') ;
                        k=k+1;
                        if max(stat)>0
                            poly1= subtract(poly1,Floe(j).poly);
                        end
                    end
                end
            end
            poly1 = subtract(poly1,polyout);
            if area(poly1)>0
                clear c0;
                c0(:,1) = poly1.Vertices(:,1); c0(:,2) = poly1.Vertices(:,2);
                c0 = c0';
                breaks = isnan(poly1.Vertices(:,1));
                I = find(breaks == 1);
                [m,n] = size(I);
                if m  == 1
                    I = [0 I length(breaks)+1];
                elseif n == 1
                    I = [0 I' length(breaks)+1];
                end
                for jj = 1:poly1.NumRegions
                    Anew = area(poly1,jj)+Anew;
                    polynew = polyshape(c0(1,I(jj)+1:I(jj+1)-1),c0(2,I(jj)+1:I(jj+1)-1));
                    cnow = (Area+Anew)/((max(c2_boundary(2,:))-min(c2_boundary(2,:)))*(max(c2_boundary(1,:))-min(c2_boundary(1,:))));
                    [Floe2(ind)] = initialize_floe_values(polynew);
                    Floe2(ind).potentialInteractions = FloeNEW(ii).potentialInteractions;
                    ind = ind+1;
                end
            end
        end
    end
    
    if ii >= length(C)
%     if ii >= length(FloeV)
        cnow = 1;
    end
    ii = ii+1;
end
warning('on',id)

Floe2=rmfield(Floe2,'potentialInteractions');
end