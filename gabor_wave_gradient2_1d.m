%This function computes the gradient of X_tilda and update it to make
%wavelet coefficients more close to Sx(wavelet coefficient of x)
%Compute gradient in frequency
function [f, g] = gabor_wave_gradient2_1d(X_tilda, Sx, psi, psi_hat)

N = length(X_tilda);
nS = size(psi_hat,2);
Sx_tilda = scattering_gabor_infreq_1d(X_tilda, psi_hat);
f = sum((Sx_tilda - Sx).^2);
g = zeros(N,1);
X_tilda_hat = fft(fftshift(X_tilda));
eps = 1e-6;
for i = 1:nS
    temp1 = fftshift(ifft(X_tilda_hat .* psi_hat(:,i)));
    for l = 1:N
        G(:,l) = circshift(psi(:,i), l-1);
    end
%     temp2 = ((real(temp1) ./ (abs(temp1) + eps)' * real(G) + (imag(temp1)./ ...
%         (abs(temp1) + eps)))' * imag(G));
    temp2 = ((real(temp1)./ abs(temp1))' * real(G) + (imag(temp1)./ ...
        abs(temp1))' * imag(G));
    g = g + 2 * (Sx_tilda(i) - Sx(i)) * temp2';
end
g = g + 2 * (Sx_tilda(end) - Sx(end)) * sign(X_tilda);
