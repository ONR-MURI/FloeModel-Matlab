function psi_hat = CalcStreamFcn(q_hat,p)
% Calculate streamfunctions of PV (q_hat) and struct containing
%q=p;
%p = real(q);

persistent Laplacian InvBT InvBC
if isempty(Laplacian)
    k = [0:p.N/2 -p.N/2+1:-1]'/p.hs;
    dX = 1i*repmat(k',[p.N 1 2]);
    dY = 1i*repmat(k,[1 p.N 2]);
    Laplacian = dX(:,:,1).^2+dY(:,:,1).^2;
    InvBT = 1./Laplacian; InvBT(1,1) = 0;
    InvBC = 1./(Laplacian-p.kd);InvBC(1,1) = 0;
end

% Invert for psi
q_bt = .5*(q_hat(:,:,1) + q_hat(:,:,2)); % q_bt = q_psi in Qi & Majda 2016 paper, bt barotropic
q_bc = .5*(q_hat(:,:,1) - q_hat(:,:,2)); % q_bc = q_tau in Qi & Majda 2016 paper, bc baroclinic
psi_bt = InvBT.*q_bt;  %psi in Qi & Majda 2016 paper
psi_bc = InvBC.*q_bc;  %tau in Qi & Majda 2016 paper
psi_hat(:,:,2) = psi_bt-psi_bc; % psi_2
psi_hat(:,:,1) = psi_bt+psi_bc; % psi_1

end