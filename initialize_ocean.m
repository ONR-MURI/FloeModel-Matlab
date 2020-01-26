function ocean=initialize_ocean()

% defining ocean currents

% ocean grid;
dXo=2000; % in meters

Xo=-70e3:dXo:70e3; Yo=-70e3:dXo:70e3; 
[Xocn, Yocn]=meshgrid(Xo,Yo);

%defining ocean streamfunction with some eddies
transport=1e3; % horizontal transport, in m^2/s (controls ocean currents) 
psi_ocean=transport*sin(2*pi*Xocn/40e3).*cos(2*pi*Yocn/50e3); 


%calculating ocean velocity field 
Uocn=zeros(size(Xocn)); Vocn=zeros(size(Xocn));
Uocn(2:end,:)=-(psi_ocean(2:end,:)-psi_ocean(1:end-1,:))/dXo; 
Vocn(:,2:end)=(psi_ocean(:,2:end)-psi_ocean(:,1:end-1))/dXo;
%adding pure divergence 
%Uocn=Uocn-0.1*Xocn/3e4;
%Vocn=Vocn-0.1*Yocn/3e4;

ocean.Xo=Xo;
ocean.Yo=Yo;
ocean.Uocn=Uocn;
ocean.Vocn=Vocn;

end