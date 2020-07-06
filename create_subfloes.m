function [subfloes,X,Y,boundingbox] = create_subfloes(floe,N,EXISTING)
%%This function takes in a polyshape and uses voronoi tesselation to break
%%it up into smaller subfloes

%Identify if there are existing points for the voronoi program or if new
%ones need to be created
if EXISTING
    X = floe.vorX; Y = floe.vorY;
    boundingbox = floe.vorbox;
else
    X = floe.Xi+floe.rmax*(2*rand(N,1)-1);
    Y = floe.Yi+floe.rmax*(2*rand(N,1)-1);
    boundingbox=[-1 ,-1; 1,-1; 1,1; -1 ,1]*sqrt(2)*floe.rmax+[floe.Xi floe.Yi];
end

%Create the new polyshapes over an entire box
worked = 1;
while worked > 0.5
    [~, b,~,~,worked] = polybnd_voronoi([X Y],boundingbox);
    k = 1;
    if worked == 1
        X = floe.Xi+floe.rmax*(2*rand(N,1)-1);
        Y = floe.Yi+floe.rmax*(2*rand(N,1)-1);
    end
end

%Find where these shapes intersect with the floe to create the new subfloes
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