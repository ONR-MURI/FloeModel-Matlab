function RHS = RHS_Spectral(q_hat,p)
% Function takes Fourier coefficients of PV (q_hat) and struct containing
% parameters (p) and evaluates RHS of 2-layer QG equations except for
% high-k dissipation. Returns Fourier coefficients of RHS.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent dX dY DX DY Laplacian InvBT InvBC
if isempty(DX)
    k = [0:p.N/2 -p.N/2+1:-1]'/p.hs;  % modes scaled to [-1e6, 1e6]
    dX = 1i*repmat(k',[p.N 1 2]);
    dY = 1i*repmat(k,[1 p.N 2]);
    Laplacian = dX(:,:,1).^2+dY(:,:,1).^2;
    InvBT = 1./Laplacian; InvBT(1,1) = 0;
    InvBC = 1./(Laplacian-p.kd);InvBC(1,1) = 0;
    k = [0:p.N/2-1 0 -p.N/2+1:-1]'/p.hs;
    dX = 1i*repmat(k',[p.N 1 2]);
    dY = 1i*repmat(k,[1 p.N 2]);
    % For the dealiased jacobian:
    k = [0:.75*p.N-1 0 -.75*p.N+1:-1]'/p.hs;
    DX = 1i*repmat(k',[1.5*p.N 1 2]);
    DY = 1i*repmat(k,[1 1.5*p.N 2]);
    clear k
end

% Invert for psi
q_bt = .5*(q_hat(:,:,1) + q_hat(:,:,2)); % q_hat is the q1 q2, q_bt = q_psi in Qi & Majda 2016 paper, bt barotropic
q_bc = .5*(q_hat(:,:,1) - q_hat(:,:,2)); % q_bc = q_tau in Qi & Majda 2016 paper, bc baroclinic
psi_bt = InvBT.*q_bt;  %psi in Qi & Majda 2016 paper
psi_bc = InvBC.*q_bc;  %tau in Qi & Majda 2016 paper
psi_hat(:,:,2) = psi_bt-psi_bc; % psi_2
psi_hat(:,:,1) = psi_bt+psi_bc; % psi_1

% calculate Ekman plus beta plus mean shear
RHS = zeros([p.N p.N 2]);
RHS(:,:,1) =-p.U*dX(:,:,1).*q_hat(:,:,1)-(p.kb + p.U*p.kd)*dX(:,:,1).*psi_hat(:,:,1);
RHS(:,:,2) = p.U*dX(:,:,1).*q_hat(:,:,2)-(p.kb - p.U*p.kd)*dX(:,:,1).*psi_hat(:,:,2) - p.r*Laplacian.*psi_hat(:,:,2);
    
% For using a 3/2-rule idealiased jacobian:
    % physical space, 3/2 grid; factor of (9/4) scales fft
    Psi_hat = zeros([1.5*p.N 1.5*p.N 2]);
    Psi_hat(1:p.N/2+1,1:p.N/2+1,:) = (9/4)*psi_hat(1:p.N/2+1,1:p.N/2+1,:);
    Psi_hat(1:p.N/2+1,p.N+2:1.5*p.N,:) = (9/4)*psi_hat(1:p.N/2+1,p.N/2+2:p.N,:);
    Psi_hat(p.N+2:1.5*p.N,1:p.N/2+1,:) = (9/4)*psi_hat(p.N/2+2:p.N,1:p.N/2+1,:);
    Psi_hat(p.N+2:1.5*p.N,p.N+2:1.5*p.N,:) = (9/4)*psi_hat(p.N/2+2:p.N,p.N/2+2:p.N,:);
    Q_hat = zeros([1.5*p.N 1.5*p.N 2]);
    Q_hat(1:p.N/2+1,1:p.N/2+1,:) = (9/4)*q_hat(1:p.N/2+1,1:p.N/2+1,:);
    Q_hat(1:p.N/2+1,p.N+2:1.5*p.N,:) = (9/4)*q_hat(1:p.N/2+1,p.N/2+2:p.N,:);
    Q_hat(p.N+2:1.5*p.N,1:p.N/2+1,:) = (9/4)*q_hat(p.N/2+2:p.N,1:p.N/2+1,:);
    Q_hat(p.N+2:1.5*p.N,p.N+2:1.5*p.N,:) = (9/4)*q_hat(p.N/2+2:p.N,p.N/2+2:p.N,:);
    % calculate u.gradq on 3/2 grid
    u = real(ifft2(-DY.*Psi_hat));
    v = real(ifft2( DX.*Psi_hat));
    qx= real(ifft2( DX.*Q_hat));
    qy= real(ifft2( DY.*Q_hat));
    jaco_real = u.*qx+v.*qy;
    % fft, 3/2 grid; factor of (4/9) scales fft
    Jaco_hat = (4/9)*fft2(jaco_real);
    % reduce to normal grid
    jaco_hat = zeros([p.N p.N 2]);
    jaco_hat(1:p.N/2+1,1:p.N/2+1,:) = Jaco_hat(1:p.N/2+1,1:p.N/2+1,:);
    jaco_hat(1:p.N/2+1,p.N/2+2:p.N,:) = Jaco_hat(1:p.N/2+1,p.N+2:1.5*p.N,:);
    jaco_hat(p.N/2+2:p.N,1:p.N/2+1,:) = Jaco_hat(p.N+2:1.5*p.N,1:p.N/2+1,:);
    jaco_hat(p.N/2+2:p.N,p.N/2+2:p.N,:) = Jaco_hat(p.N+2:1.5*p.N,p.N+2:1.5*p.N,:);
    
% Put it all together
RHS = RHS - jaco_hat;

%% Two-layer QG above; now we add the effect of heat flux from the ocean/sea ice
% for the lower layer (q_2) in this code is -1/taub*(\psi_2 -
% \psi_1)*kd - Tsrf/dt); 
HF(:,:,1) = -p.taub*( (psi_hat(:,:,2)-psi_hat(:,:,1))*p.kd + p.Tsf*p.Tc/p.dz );
HF(:,:,2) = -p.taub*( (psi_hat(:,:,1)-psi_hat(:,:,2))*p.kd - p.Tsf*p.Tc/p.dz );

RHS = RHS + 0*HF;

end