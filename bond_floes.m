function Floe=bond_floes(Floe)

x=cat(1,Floe.Xi);
y=cat(1,Floe.Yi);
rmax=cat(1,Floe.rmax);
alive=cat(1,Floe.alive);

for i=1:length(Floe)
    Floe(i).bonds=[];
    for j=1:length(Floe)
        
        if j>i && alive(j) && alive(i) && ((x(i)-x(j))^2 + (y(i)-y(j))^2)<(rmax(i)+rmax(j))^2 % if floes are potentially overlapping
        
        intersection=intersect(Floe(i).poly,Floe(j).poly);
        
        if ~isempty (intersection)
            
            k=length(Floe(i).bonds)+1;
            
            Floe(i).bonds(k).FloeNum=j;
            
            Floe(i).bonds(k).Points=[intersection.Vertices intersection.Vertices];
            
            Floe(i).bonds(k).Strength=1; % this can hold memory of stains

        
            k=length(Floe(j).bonds)+1;
            
            Floe(j).bonds(k).FloeNum=i;
            
            Floe(j).bonds(k).Points=[intersection.Vertices intersection.Vertices];
            
            Floe(j).bonds(k).Strength=1; % this can hold memory of strains  
        
        end
        
        end 
        
    end
    
end

end