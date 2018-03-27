function Sx = scattering_gabor_infreq_1d(X, psi_hat)

N = length(X);
ncol = size(psi_hat,2);
%Sx = zeros(ncol, 1);
Sx = zeros(ncol + 1, 1);
X_hat = fft(X);

for ind = 1:ncol
    temp = ifft(X_hat .* psi_hat(:,ind));
    Sx(ind) = sum(abs(temp));
end

Sx(ncol + 1) = sum(abs(X));


