% Linear Gassian model for atmosphere winds fields

clear all; close all; clc;

load('AtmThetaN100.mat','theta');
load('AtmOmegaN100.mat','omega');
load('AtmMuN100.mat','mu');
load('AtmSigmaN100.mat','sigma');
load('AtmPsi0N100.mat','psi0');
load('AtmPsiNfftvarN100dt2kNt10k.mat','psiNfftvar'); 
ufftvar = psiNfftvar;

Nl = 5; N =100; h = 1/N; x = 0:h:1; dt = 0.01; Nt = 10000; t = 0:dt:Nt*dt;
if isempty(dir('dataAtmOU')); mkdir('dataAtmOU'); end

fig = figure; figg = figure;
uf = zeros(N,N,Nt); uf(:,:,1) = psi0;
for jt=2:1001 % or longer with Nt; Can be done by parfor    
    z = 1*jt*dt;
    
    for j=1:N
        for k=1:N
            if((j<Nl+1 || j>N-Nl) && (k<Nl+1 || k>N-Nl))
                uf(j,k,jt) = uf(j,k,jt-1) - (theta(j,k) + 1i*omega(j,k)).*uf(j,k,jt-1).*dt...
                    + mu(j,k).*dt + sigma(j,k).*(randn(1) + 1i*randn(1)).*sqrt(dt/2);
            end
        end
    end

    % only use noise for small scale modes
    for j=1:N
        for k=1:N
            if(~((j<Nl+1 || j>N-Nl) && (k<Nl+1 || k>N-Nl)))
                uf(j,k,jt) = mu(j,k) + sqrt(ufftvar(j,k)).*(randn(1) + 1i*randn(1)).*sqrt(dt/2);
            end
        end
    end
    
    uf(1,1,jt) = mu(1,1) + sqrt(ufftvar(1,1)).*(randn(1) + 1i*randn(1)).*sqrt(dt/2);
    
    if(mod(jt,10)==0)
    sol = real(ifft2(uf(:,:,jt)));
    figure(fig);
    colormap default
    imagesc([0,1], [0,1], sol,'CDataMapping','scaled');
    %imagesc([0,1], [0,1], real(ifft2(ufft(:,:,10*jt))),'CDataMapping','scaled'); % synthetic sols
    caxis manual
    caxis([-4e6 4e6]);
    title(['Time = ' num2str((jt/10-1)*dt*1e5) ' s']);
    colorbar
    set(gca,'Ydir','normal');
    drawnow;
    saveas(fig,['./dataAtmOU/' num2str(jt/10,'%03.f')],'jpg');
    
    winds = UpdateWinds(sol,1e6,1e6,2e4);
    figure(figg);
    dn=4; % plot every dn'th velocity vector
    quiver(winds.X(1:dn:end),winds.Y(1:dn:end),winds.U(1:dn:end,1:dn:end),winds.V(1:dn:end,1:dn:end));
    
    title(['Time = ' num2str((jt/10-1)*dt*1e5) ' s']);
    %colormap summar
    %colormap('gray'); caxis([0 1]);
    axis([-1e6 1e6 -1e6 1e6])
    xlabel('m');ylabel('m');
    set(gca,'Ydir','normal');
    saveas(figg,['./dataAtmOU/' num2str(jt/10,'w%03.f') '.jpg'],'jpg');
    end
end