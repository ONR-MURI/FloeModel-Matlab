function [floes] = FloeSimplify(floe, ddx,SUBFLOES)
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

%Check for holes in polyshapes
if floe.poly.NumHoles > 0
    if abs(area(floe.poly,2)) > 1e7
        xx = 1;
        xx(1) = [1 2];
    end
    floe.poly = rmholes(floe.poly);
end

%Set function to define the tolerance for Douglas-Peuker algorithm
SimpMin = @(A) 0.5*log10(A)^2.5;

%Create new simplified polyshape
floenew = floe;
vertnew = DouglasPeucker(floe.poly.Vertices,SimpMin(floe.area));
if length(vertnew) >= length(floe.poly.Vertices)
    vertnew = DouglasPeucker(floe.poly.Vertices,2*SimpMin(floe.area));
end
polynew = polyshape(vertnew);
if sum(area(intersect(floe.poly,polynew)))/floe.area < 0.5 && area(polynew) > 1e6
    xx = 1;
    xx(1) = [1 2];
end

%Align center of old polygon with the enw one
[x1,y1] = centroid(polynew);
dx = floe.Xi-x1;
dy = floe.Yi-y1;
if area(intersect(floe.poly,polynew))/floe.area > 0.95
    polynew = translate(polynew,[dx, dy]);
end

%Check if simplifcation led to polygon having multiple regions
polyout = sortregions(polynew,'area','descend');
R = regions(polyout);
R = R(area(R)>1e4);

%% Calculate the new properties associated with this floe since it has a new shape
if ~isempty(R)
    for ii = 1:length(R)
        poly1new = R(ii);
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
        k = 1;
        
        clear subfloes
        if SUBFLOES
            [~, b,~,~,worked] = polybnd_voronoi([X Y],boundingbox);
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
        
        if floenew.poly.NumRegions > 1
            xx = 1;
            xx(1) = [1 2];
        end
        
        %populate new floes with heights of exisiting floes
        N = length(subfloes);
        floenew.SubFloes = [];
        areaS = zeros(N,1);
        inertia = zeros(N,1);
        centers = zeros(N,2);
        [floex, floey] = centroid([floe.SubFloes.poly]);
        for kk = 1:N
            [Xi,Yi] = centroid(subfloes(kk));
            [I,~] = dsearchn([floex',floey'],[Xi,Yi]);
            centers(kk,:) = [Xi,Yi];
            floenew.SubFloes(kk).poly = subfloes(kk);
            floenew.SubFloes(kk).h = floe.SubFloes(I).h;
            areaS(kk) = area(floenew.SubFloes(kk).poly);
            if floenew.SubFloes(kk).poly.NumHoles > 0
                breaks = isnan(floenew.SubFloes(kk).poly.Vertices(:,1));
                I = find(breaks == 1);
                I = [0 I' length(breaks)+1];
                inertia2 = zeros(floenew.SubFloes(kk).poly.NumHoles,1);
                centers2 = zeros(floenew.SubFloes(kk).poly.NumHoles,2);
                for jj = 1:length(I)-1 
                    [Xi,Yi] = centroid(polyshape(floenew.SubFloes(kk).poly.Vertices(I(jj)+1:I(kk+1)-1,:)));
                    centers2(kk,:) = [Xi,Yi];
                    inertia2(kk) = PolygonMoments(floenew.SubFloes(kk).poly.Vertices(I(jj)+1:I(kk+1)-1,:),floenew.SubFloes(kk).h);
                end
                inertia(kk) = sum(inertia2'+floenew.SubFloes(kk).h.*sqrt((centers2(:,1)-centers2(1,1)).^2+(centers(:,2)-centers2(1,2)).^2));
            else
                inertia(kk) = PolygonMoments(floenew.SubFloes(kk).poly.Vertices,floenew.SubFloes(kk).h);
            end
        end
        
        %%
        floenew.mass = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h));
        floenew.Xm = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h).*centers(:,1))./floenew.mass;
        floenew.Ym = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h).*centers(:,2))./floenew.mass;
        floenew.inertia_moment = sum(inertia+cat(1,floenew.SubFloes.h).*sqrt((centers(:,1)-floenew.Xm).^2+(centers(:,2)-floenew.Ym).^2));
        floenew.rmax = sqrt(max(sum((floenew.poly.Vertices' - [floenew.Xi; floenew.Yi]).^2,1)));        

        floes = [floes floenew];
    end
else
    floes = floe;
    floes.alive = 0;
end

warning('on',id)
warning('on',id2)
warning('on',id3)
end

