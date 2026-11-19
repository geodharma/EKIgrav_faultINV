%% Generate synthetic data (additive Gaussian noise, sigma = 5% of anomaly range)
clear; close all; clc;
rng(20);                                   % reproducible noise realization
x = (0:5:300)';  N = numel(x);
m_true1 = [150; 15;  9; 25; 180];
m_true2 = [100; 15; 45; 180];
m_true3 = [100; 16; 90];
noise_level = 0.05;                        % 5% of the anomaly range

g_true = [sesar_Eq1(m_true1, x')', sesar_Eq2(m_true2, x')', sesar_Eq3(m_true3, x')'];
g_noisy = zeros(N, 3);
sigma_n = zeros(1, 3);
for k = 1:3
    sigma_n(k)    = noise_level * (max(g_true(:,k)) - min(g_true(:,k)));
    g_noisy(:,k)  = g_true(:,k) + sigma_n(k) * randn(N, 1);
    noise_k       = g_noisy(:,k) - g_true(:,k);
    fprintf('Model %d: sigma_n = %.3f mGal | RMS noise = %.3f mGal (%.1f%% of range) | SNR = %.2f dB\n', ...
            k, sigma_n(k), rms(noise_k), 100*rms(noise_k)/(max(g_true(:,k))-min(g_true(:,k))), ...
            20*log10(std(g_true(:,k))/std(noise_k)));
end

% Quick-look plots
figure('Color','w');
for k = 1:3
    subplot(3,1,k); plot(x, g_true(:,k), 'k-', 'LineWidth', 1.5); hold on;
    plot(x, g_noisy(:,k), 'r.', 'MarkerSize', 12);
    ylabel('\Deltag (mGal)'); title(sprintf('Model %d', k));
end
xlabel('Distance (m)');

% Save (overwrite old files)
%writematrix([x g_true],  'data_fault_cleanNEW.csv');
%writematrix([x g_noisy], 'data_fault_noisyNEW.csv');
%writematrix(sigma_n,     'data_fault_sigma.csv');   % sigma_n per model, for sd_kf in EKI