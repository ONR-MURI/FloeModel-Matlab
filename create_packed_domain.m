function Floe= create_packed_domain(Floe,dL)

%dL=4e4;
%Floe=load('Floe_packed.mat','Floe'); Floe=Floe.Floe;
Floe1=Floe;
Floe2=Floe;
Floe3=Floe; 
Floe4=Floe;

for i=1:length(Floe)
    
    Floe1(i).Xi = Floe(i).Xi+dL;
    Floe1(i).Yi = Floe(i).Yi+dL;
    Floe1(i).poly.Vertices= Floe(i).poly.Vertices + [dL dL];


    Floe2(i).Xi = Floe(i).Xi-dL;
    Floe2(i).Yi = Floe(i).Yi+dL;
    Floe2(i).poly.Vertices= Floe(i).poly.Vertices + [-dL dL];


    Floe3(i).Xi = Floe(i).Xi+dL;
    Floe3(i).Yi = Floe(i).Yi-dL;
    Floe3(i).poly.Vertices= Floe(i).poly.Vertices + [dL -dL];
    
    Floe4(i).Xi = Floe(i).Xi-dL;
    Floe4(i).Yi = Floe(i).Yi-dL;
    Floe4(i).poly.Vertices= Floe(i).poly.Vertices + [-dL -dL];

end

Floe=[Floe1 Floe2 Floe3 Floe4];

end
