function Floe= create_MIZdomain()

load('PackedFloesFullDomain.mat','Floe'); 

Floe1=Floe;
Floe2=Floe;
Floe3=Floe;
Floe4=Floe;
dL=130e3;
for i=1:length(Floe)
    
    Floe1(i).Xi = Floe(i).Xi-1.5*dL;
    Floe1(i).Yi = Floe(i).Yi+dL/2;
    
    Floe2(i).Xi = Floe(i).Xi-0.5*dL;
    Floe2(i).Yi = Floe(i).Yi+dL/2;

    Floe3(i).Xi = Floe(i).Xi+0.5*dL;
    Floe3(i).Yi = Floe(i).Yi+dL/2;

    Floe4(i).Xi = Floe(i).Xi+1.5*dL;
    Floe4(i).Yi = Floe(i).Yi+dL/2;
    
end

Floe=[Floe1 Floe2 Floe3 Floe4];

end
