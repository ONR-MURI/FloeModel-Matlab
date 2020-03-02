
function Floes=fracture_floe(Floe,N)


Floes=[]; k=1;

for p=1:length(Floe)
    
    floe=Floe(p);
    Anew = 0;
    
    while Anew/area(floe.poly) < 0.9
        k2 = k;
        X = floe.Xi+floe.rmax*(2*rand(N,1)-1);
        Y = floe.Yi+floe.rmax*(2*rand(N,1)-1);
        
        
        boundingbox=[-1 ,-1; 1,-1; 1,1; -1 ,1]*floe.rmax+[floe.Xi floe.Yi];
        [~, b] = polybnd_voronoi([X Y],boundingbox);
        
        
        for i=1:length(b)
            a=intersect(floe.poly,polyshape(b{i}));
            if a.NumRegions==1,
                if k2==1 Floes=initialize_floe_values(a); k2=k2+1;
                else Floes(k2)=initialize_floe_values(a);k2=k2+1;
                end
            end
            
        end
        
        Anew = sum(area([Floes.poly]));
        
    end
    k = k2;
    
end

% figure; 
% plot(floe.poly)
% hold on;
% plot(X,Y,'x');
% plot([Floes.poly])



end