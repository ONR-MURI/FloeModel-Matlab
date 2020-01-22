function Floe= create_polygons(Floe)

for i=1:length(Floe)            
  
    Floe(i).poly=polyshape(Floe(i).c_alpha(1,:)+Floe(i).Xi,Floe(i).c_alpha(2,:)+Floe(i).Yi); 
    
    Floe(i).r_max=sqrt(max(sum((Floe(i).poly.Vertices' - [Floe(i).Xi; Floe(i).Yi]).^2,1))); % exact Rmax.
    
%    Floe(i).c0 = Floe(i).poly.Vertices';
%    Floe(i).c_alpha =Floe(i).c0;
    
end

end