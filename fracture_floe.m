
function Floes=fracture_floe(Floe,N)

Floes=[]; k=1;

for p=1:length(Floe)
    
    floe=Floe(p);


X = floe.Xi+floe.rmax*(2*rand(N,1)-1);
Y = floe.Yi+floe.rmax*(2*rand(N,1)-1);


boundingbox=[-1 ,-1; 1,-1; 1,1; -1 ,1]*floe.rmax+[floe.Xi floe.Yi];
[~, b] = polybnd_voronoi([X Y],boundingbox);


for i=1:length(b)
    a=intersect(floe.poly,polyshape(b{i})); 
    if a.NumRegions==1, 
        if k==1 Floes=initialize_floe_values(a); k=k+1;
        else Floes(k)=initialize_floe_values(a);k=k+1; 
        end
    end
    
end




end

% figure; 
% plot(floe.poly)
% hold on;
% plot(X,Y,'x');
% plot([Floes.poly])



end