function [X_bar, SX_bar, Sx] = gaborwave_synth_greedy_1d(X, S, xi, x0)
% This function synthesizes one dimensional vector by minimizing Gabor
% difference using greedy algorithm
%
% Needed functions:
%   scattering_gabor_inspace_1d.m, gabor_wave_diff_1d.m 
%
% Args:
%   X: the vector wanted to synthesize
%   S: vector of scales
%   xi: frequency parameter
%   k: number of iterations
% Returns:
%   X_bar: the synthetic vector which is similar to X
%   SX_bar: scattering features of synthesized signal
%   SX: scattering features of target signal
%   This function also plots the original vector and produced vector to
%   make comparison.

N = size(X, 1);

u = [1:N];
g = zeros(N, length(S));
g_hat = zeros(N, length(S));
%g_bar = zeros(N, length(S));
for ind = 1:length(S)
    s = S(ind);
%     g(1:s, ind) = 1/s .* exp(i * xi * u(1:s))';
    g(1:s, ind) = exp(i * xi * u(1:s))';
   % g_bar((N-s+1):N, ind) = - exp(i * xi * u((N-s+1):N))'; 
        %negative because we need (-1)^k where xi = k/N * pi, here k = 41
    g_hat(:,ind) = fft(g(:,ind));
   % g_bar_hat(:,ind) = fft(g_bar(:,ind));
end

%Sx = scattering_gabor_inspace_1d(X, g);
Sx = scattering_gabor_infreq_1d(X, g_hat);
%define the function we need to minimize sum squares.

%initialize X0
if nargin < 4
    X0 = randn(N,1);
else
    X0 = x0;
end

options = optimoptions('fminunc',...%'Display','iter',...
    'MaxFunctionEvaluations',10000000,'SpecifyObjectiveGradient',true, ...
    'MaxIterations', 50000, 'FunctionTolerance', 1e-8,'StepTolerance', 1e-8);


SX_bar = cell(length(S), 1);
X_bar = zeros(N,length(S));
fun = @(X_tilda)gabor_wave_diff_1d(X_tilda, Sx, g_hat);

for ind = 1:length(S)
    %fun = @(X_tilda)gabor_wave_diff_1d(X_tilda, [Sx(1:ind);Sx(end)], g_hat(:,1:ind));
    fun = @(X_tilda)gabor_wave_gradient_infreq_1d(X_tilda, [Sx(1:ind);Sx(end)], ...
        g(:,1:ind), g_hat(:,1:ind));
    SX_bar{ind} = scattering_gabor_infreq_1d(X0, g_hat(:,1:ind));
    c = 0;
    C = 0;
    while norm(SX_bar{ind} - [Sx(1:ind);Sx(end)])^2 > 1e-4
        X_bar(:,ind) = fminunc(fun, X0, options);
        SX_bar{ind} = scattering_gabor_infreq_1d(X_bar(:,ind), g_hat(:,1:ind));
        if c < 10
            X0 = X_bar(:,ind) + randn(N,1)/50;
            c = c + 1;
        elseif C <= 1
            X0 = X_bar(:,max(ind-1,1)) + randn(N,1)/5;
            c = 0;
            C = C + 1;
        else
            break;
        end
    end
    X0 = X_bar(:,ind);
end
    
