function [floes] = FloeSimplify(floe, ddx,SUBFLOES)
%Take polyshape with a lot of vertices and simplify it to have fewer
%vertices
%% Remap the main polygon to a shape with fewer vertices
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id2 = 'MATLAB:polyshape:boolOperationFailed';
warning('off',id2)
floes = [];
rho_ice = 920;
if floe.poly.NumHoles > 0
    if abs(area(floe.poly,2)) > 1e7
        xx = 1;
        xx(1) = [1 2];
    end
    floe.poly = rmholes(floe.poly);
end
SimpMin = @(A) 0.5*log10(A)^2.5;

floenew = floe;
n=(fix(floe.rmax/ddx)+1); n=ddx*(-n:n);
[Xg, Yg]= meshgrid(n+floe.Xi, n+floe.Yi);
P = [Xg(:) Yg(:)];

 vertnew = DouglasPeucker(floe.poly.Vertices,SimpMin(floe.area));
 polynew = polyshape(vertnew);
 if sum(area(intersect(floe.poly,polynew)))/floe.area < 0.5 && area(polynew) > 3500
     xx = 1;
     xx(1) = [1 2];
 end
% [k,~] = dsearchn(P,floe.poly.Vertices);
% polynew = polyshape(P(k,1),P(k,2));
[x1,y1] = centroid(polynew);
dx = floe.Xi-x1;
dy = floe.Yi-y1;
if area(intersect(floe.poly,polynew))/floe.area > 0.95
    polynew = translate(polynew,[dx, dy]);
% elseif area(intersect(floe.poly,polynew))/floe.area < 0.75 && floe.area > 1e7
%     xx = 1;
%     xx(1) = [1 2];
end
polyout = sortregions(polynew,'area','descend');
R = regions(polyout);
R = R(area(R)>1e4);
%% 

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
%         worked = 1;
%         while worked > 0.5
            
%             if worked == 1
%                 X = floe.Xi+floe.rmax*(2*rand(N,1)-1);
%                 Y = floe.Yi+floe.rmax*(2*rand(N,1)-1);
%             end
%         end
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
        for ii = 1:N
            [Xi,Yi] = centroid(subfloes(ii));
            [I,~] = dsearchn([floex',floey'],[Xi,Yi]);
            centers(ii,:) = [Xi,Yi];
            %     poly = intersect(subfloes(ii),[floe.SubFloes.poly]);
            %     [~,I] = min(abs(area(poly)/area(subfloes(ii))-1));
            floenew.SubFloes(ii).poly = subfloes(ii);
            floenew.SubFloes(ii).h = floe.SubFloes(I).h;
            areaS(ii) = area(floenew.SubFloes(ii).poly);
            if floenew.SubFloes(ii).poly.NumHoles > 0
                breaks = isnan(floenew.SubFloes(ii).poly.Vertices(:,1));
                I = find(breaks == 1);
                I = [0 I' length(breaks)+1];
                inertia2 = zeros(floenew.SubFloes(ii).poly.NumHoles,1);
                centers2 = zeros(floenew.SubFloes(ii).poly.NumHoles,2);
                for kk = 1:length(I)-1 
                    [Xi,Yi] = centroid(polyshape(floenew.SubFloes(ii).poly.Vertices(I(kk)+1:I(kk+1)-1,:)));
                    centers2(kk,:) = [Xi,Yi];
                    inertia2(kk) = PolygonMoments(floenew.SubFloes(ii).poly.Vertices(I(kk)+1:I(kk+1)-1,:),floenew.SubFloes(ii).h);
                end
                inertia(ii) = sum(inertia2'+floenew.SubFloes(ii).h.*sqrt((centers2(:,1)-centers2(1,1)).^2+(centers(:,2)-centers2(1,2)).^2));
%                 xx = 1;
%                 xx(1) = [1 2];
            else
                inertia(ii) = PolygonMoments(floenew.SubFloes(ii).poly.Vertices,floenew.SubFloes(ii).h);
            end
        end
        
        %%
        floenew.mass = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h));
        floenew.Xm = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h).*centers(:,1))./floenew.mass;
        floenew.Ym = sum(rho_ice*areaS.*cat(1,floenew.SubFloes.h).*centers(:,2))./floenew.mass;
        floenew.inertia_moment = sum(inertia+cat(1,floenew.SubFloes.h).*sqrt((centers(:,1)-floenew.Xm).^2+(centers(:,2)-floenew.Ym).^2));
        floenew.rmax = sqrt(max(sum((floenew.poly.Vertices' - [floenew.Xi; floenew.Yi]).^2,1)));
        
        % if sum(floenew.h.*areaS)/floenew.area<sum([floe.SubFloes.h].*area([floe.SubFloes.poly]))/floe.area && length(floenew.SubFloes)<= length(floe.SubFloes)
        %     xx = 1;
        %     xx(1) = [1 2];
        % end
        floes = [floes floenew];
    end
else
    floes = floe;
    floes.alive = 0;
end

warning('on',id)
warning('on',id2)
end

