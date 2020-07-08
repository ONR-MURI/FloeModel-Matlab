function winds = InitializeWinds(transport, Lx, Ly,dL)
% Initializing winds velocity field 

X = -Lx:dL:Lx-dL; Y = -Ly:dL:Ly-dL;
[Xatm, Yatm] = meshgrid(X,Y); 
psiW = transport*sin(2*pi*Xatm/Lx).*cos(1*pi*Yatm/Ly); 

U = zeros(size(X)); V = zeros(size(Y));
U(2:end,:) = -(psiW(2:end,:)-psiW(1:end-1,:))/dX; 
V(:,2:end) = (psiW(:,2:end)-psiW(:,1:end-1))/dX;

winds.X = X; winds.Y = Y;
winds.U = U; winds.V = V;

end