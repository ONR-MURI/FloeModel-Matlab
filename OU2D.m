function [tpsi] = OU2D(theta, omega, sigma, mu, psivar, psi, N, dt)

tpsi = zeros(N); Nl = 5; % Taking five dominated modes
dt = dt*1e-5; % Time scaling from OU modeling

for j=1:N
    for k=1:N
        if((j<Nl+1 || j>N-Nl) && (k<Nl+1 || k>N-Nl))
            tpsi(j,k) = psi(j,k) - (theta(j,k) + 1i*omega(j,k)).*psi(j,k).*dt...
                + mu(j,k).*dt + sigma(j,k).*(randn(1) + 1i*randn(1)).*sqrt(dt/2);
        end
    end
end

for j=1:N
    for k=1:N
        if(~((j<Nl+1 || j>N-Nl) && (k<Nl+1 || k>N-Nl)))
            tpsi(j,k) = mu(j,k) + sqrt(psivar(j,k)).*(randn(1) + 1i*randn(1)).*sqrt(dt/2);
        end
    end
end

tpsi(1,1) = mu(1,1) + sqrt(psivar(1,1)).*(randn(1) + 1i*randn(1)).*sqrt(dt/2);
end