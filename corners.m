function [Floe] = corners(Floe,Nb)
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
for ii = 1+Nb:length(Floe)
    if ~isempty(Floe(ii).interactions)
        inter = Floe(ii).interactions(:,1);
%         inter = inter(inter<=N0);
        Xi = Floe(ii).interactions(inter<=N0,4);%x(inter);
        Yi = Floe(ii).interactions(inter<=N0,5);%y(inter);
        d = sqrt((Xi-Floe(ii).Xi).^2-(Yi-Floe(ii).Yi).^2);
        
        poly = polyshape(Floe(ii).c_alpha');
        if poly.NumRegions == 1 && poly.NumHoles == 0
            if norm(poly.Vertices(1,:)-poly.Vertices(end,:)) == 0
                poly.Vertices(end,:) = [];
            end
            angles = polyangles(poly.Vertices(:,1),poly.Vertices(:,2));
            Anorm = 180-360/length(angles);
            keep1=(rand(length(angles),1)>angles/Anorm);
            da = zeros(length(angles),1);
            for jj = 1:length(angles)
                if ~keep1(jj)
                    d2 = sqrt((Xi-poly.Vertices(jj,1)).^2-(Yi-poly.Vertices(jj,2)).^2);
                    if max(d-d2)>0
                        da(jj) = 1;
                    end
                end
            end
            
            keep = logical(keep1 + da);
            
            if sum(keep)>0
                [ang,~] = min(angles(keep));
                I = find(angles == ang);
                if length(I) > 1
                    tick = 0; count = 1;
                    while tick < 1
                        if keep(I(count)) == 1
                            jj = I(count); tick = 2;
                        else
                            count = count + 1;
                        end
                    end
                else
                    jj = I;
                end
                KeepF(ii) = 0;
                fracturedFloes = frac_corner(Floe(ii),jj,poly);
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

