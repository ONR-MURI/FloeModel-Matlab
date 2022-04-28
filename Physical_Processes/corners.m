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
        Xi = Floe(ii).interactions(inter<=N0,4);
        Yi = Floe(ii).interactions(inter<=N0,5);
        d = sqrt((Xi-Floe(ii).Xi).^2-(Yi-Floe(ii).Yi).^2);
        
        poly = polyshape(Floe(ii).c_alpha');
        polytrue = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi, Floe(ii).Yi]);
        if poly.NumRegions == 1 && poly.NumHoles == 0
            if norm(poly.Vertices(1,:)-poly.Vertices(end,:)) == 0
                poly.Vertices(end,:) = [];
            end
            angles = polyangles(poly.Vertices(:,1),poly.Vertices(:,2));
            Anorm = 180;%180-360/length(angles);
            break1=(rand(length(angles),1)>angles/Anorm);
            da = zeros(length(angles),1);
            [break2,~] = dsearchn(polytrue.Vertices,[Xi Yi]);
            da(break2) = 1;
            
            grind = logical(break1 + da==2);
            
            if sum(grind)>0
                KeepF(ii) = 0;
                fracturedFloes = frac_corner(Floe(ii),grind,poly);
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
if isfield(floenew,'potentialInteractions')
    floenew=rmfield(floenew,{'potentialInteractions'});
end
if isfield(Floe,'potentialInteractions')
    Floe=rmfield(Floe,{'potentialInteractions'});
end

Floe = [Floe(logical(KeepF)) floenew];

warning('on',id)
warning('on',id3)

end

