function Floe=RemoveFullyOverlappingFloes(Floe)

%tic; Floe= create_polygons(Floe);toc;

x=cat(1,Floe.Xi);
y=cat(1,Floe.Yi);
rmax=cat(1,Floe.rmax);
alive=cat(1,Floe.alive);

for i=1:length(Floe)
    

    for j=1:length(Floe)
        
        if j>i && alive(j) && alive(i) && ((x(i)-x(j))^2 + (y(i)-y(j))^2)<(rmax(i)+rmax(j))^2 % if floes are potentially overlapping
                      
            a=area(intersect(Floe(i).poly,Floe(j).poly));
            
            if  a > 0.25*Floe(j).area
                Floe(j).alive=0;                
                disp(['Floe' num2str(i) ',' num2str(j) 'killed']);
            end
            
            if  a > 0.25*Floe(i).area
                Floe(i).alive=0;                
                disp(['Floe' num2str(i) ',' num2str(j) 'killed']);
            end
            
        end
    end
end

%Floe=Floe(boolean(cat(1,Floe.alive)));

end