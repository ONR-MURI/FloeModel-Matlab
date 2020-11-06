function [Floe] = corners(Floe)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

x=cat(1,Floe.Xi);
y=cat(1,Floe.Yi);
N0 = length(Floe);
floenew = [];
KeepF = ones(length(Floe),1);
for ii = 1:length(Floe)
    if ~isempty(Floe(ii).interactions)
        inter = Floe(ii).interactions(:,1);
        inter = inter(inter<=N0);
        Xi = x(inter);
        Yi = y(inter);
        d = sqrt((Xi-Floe(ii).Xi).^2-(Yi-Floe(ii).Yi).^2);
        
        poly = polyshape(Floe(ii).c_alpha');
        if poly.NumRegions>1
            polyout = sortregions(poly,'area','descend');
            R = regions(polyout);
            poly = R(1);
        end
        if poly.NumHoles==0
            angles = polyangles(poly.Vertices(:,1),poly.Vertices(:,2));
            Anorm = 180-360/length(angles);
            keep1=(rand(length(angles),1)>angles/Anorm);
            da = zeros(length(angles),1);
            for jj = 1:length(angles)
                if ~keep1(jj)
                    d2 = sqrt((Xi-poly.Vertices(jj,1)).^2-(Yi-poly.Vertices(jj,2)).^2);
                    if max(d2-d)>0
                        da(jj) = 1;
                    end
                end
            end
            
            keep = logical(keep1 + da);
            
            if sum(keep)>0
                [~,I] = min(angles(keep));
                KeepF(ii) = 0;
                fracturedFloes = frac_corner(Floe(ii),I,poly);
                floenew=[floenew fracturedFloes];
            end
        end
    end
end
if isfield(floenew,'poly')
    floenew=rmfield(floenew,{'poly'});
end
if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
end
Floe = [Floe(logical(KeepF)) floenew];

warning('on',id)
warning('on',id3)

end

