
function [subfloes,X,Y,boundingbox] = create_subfloes(floe,N,EXISTING)

if EXISTING
    X = floe.vorX; Y = floe.vorY;
    boundingbox = floe.vorbox;
else
    X = floe.Xi+floe.rmax*(2*rand(N,1)-1);
    Y = floe.Yi+floe.rmax*(2*rand(N,1)-1);
    boundingbox=[-1 ,-1; 1,-1; 1,1; -1 ,1]*sqrt(2)*floe.rmax+[floe.Xi floe.Yi];
end

[~, b] = polybnd_voronoi([X Y],boundingbox);
k = 1;

for i=1:length(b)
    a=regions(intersect(floe.poly,polyshape(b{i})));
    if ~isempty(area(a))
        for j=1:length(a)
            subfloes(k) = a(j);
            k = k+1;
        end
    end
end

end