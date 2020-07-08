function winds = UpdateWinds(psiW, Lx, Ly,dL)
% Calculate wind velocities from streamfunction psiW

X = -Lx:dL:Lx-dL; Y = -Ly:dL:Ly-dL; 
U = zeros(length(X)); V = zeros(length(Y));

U(2:end,:) = -(psiW(2:end,:)-psiW(1:end-1,:))/dL; 
V(:,2:end) = (psiW(:,2:end)-psiW(:,1:end-1))/dL;

winds.X = X; winds.Y = Y;
winds.U = U; winds.V = V;

end
