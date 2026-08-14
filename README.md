# Departmental-laboratory-project-
% Experiment 3: Sampling, Aliasing, and Sinc Reconstruction
clear; close all; clc;

%% 1. Signal & Simulation Parameters
f1 = 10;                % Frequency 1 (Hz)
f2 = 40;                % Frequency 2 (Hz)
f_max = max(f1, f2);    % Maximum frequency (40 Hz)
f_nyquist = 2 * f_max;  % Nyquist rate = 80 Hz

t_duration = 0.5;       % Signal duration in seconds

% High-resolution continuous reference signal
fs_ref = 5000;          % 5 kHz plotting resolution
t_ref = 0 : 1/fs_ref : (t_duration - 1/fs_ref);
x_ref = sin(2*pi*f1*t_ref) + 0.5*sin(2*pi*f2*t_ref);

%% 2. Define Test Sampling Frequencies
fs_cases = [200, 80, 50]; 
case_names = {
    'Oversampling (f_s > 2f_{max})', ...
    'Nyquist Rate (f_s = 2f_{max})', ...
    'Undersampling (f_s < 2f_{max})'
};

%% 3. Run Simulation & Plot
figure('Name', 'Experiment 3: Sampling and Reconstruction', 'Position', [100, 100, 1200, 800]);

for i = 1:length(fs_cases)
    fs = fs_cases(i);
    Ts = 1 / fs;
    
    % Step A: Discrete Sampling
    t_sampled = 0 : Ts : (t_duration - Ts);
    x_sampled = sin(2*pi*f1*t_sampled) + 0.5*sin(2*pi*f2*t_sampled);
    
    % Step B: Ideal Whittaker-Shannon Sinc Interpolation
    % Matrix of time differences: (N_ref x N_sampled)
    t_diff = t_ref(:) - t_sampled(:)';
    sinc_matrix = sinc(t_diff / Ts); % Note: MATLAB sinc(u) = sin(pi*u)/(pi*u)
    x_recon = (sinc_matrix * x_sampled(:))';
    
    % Step C: Reconstruction Error & RMSE
    recon_error = x_ref - x_recon;
    rmse = sqrt(mean(recon_error.^2));
    
    % Step D: Frequency Spectrum Computation (FFT)
    L = length(t_ref);
    f_axis = (0 : (L/2)) * (fs_ref / L);
    
    % FFT of Reference
    X_ref_fft = fft(x_ref);
    P2_ref = abs(X_ref_fft / L);
    P1_ref = P2_ref(1:L/2+1);
    P1_ref(2:end-1) = 2 * P1_ref(2:end-1);
    
    % FFT of Reconstructed Signal
    X_recon_fft = fft(x_recon);
    P2_recon = abs(X_recon_fft / L);
    P1_recon = P2_recon(1:L/2+1);
    P1_recon(2:end-1) = 2 * P1_recon(2:end-1);
    
    % --- Column 1: Time-Domain Waveforms ---
    subplot(3, 3, (i-1)*3 + 1);
    plot(t_ref, x_ref, 'k--', 'LineWidth', 1.2); hold on;
    stem(t_sampled, x_sampled, 'r', 'filled', 'LineWidth', 1);
    plot(t_ref, x_recon, 'b-', 'LineWidth', 1.1); hold off;
    title(sprintf('%s (f_s = %d Hz)', case_names{i}, fs));
    xlabel('Time (s)'); ylabel('Amplitude');
    xlim([0, 0.2]); % Zoomed view
    grid on;
    if i == 1
        legend('Original x(t)', 'Samples', 'Reconstructed', 'Location', 'northeast');
    end
    
    % --- Column 2: Reconstruction Error ---
    subplot(3, 3, (i-1)*3 + 2);
    plot(t_ref, recon_error, 'm-', 'LineWidth', 1.1);
    title(sprintf('Error (RMSE = %.4f)', rmse));
    xlabel('Time (s)'); ylabel('Error');
    xlim([0, 0.2]);
    grid on;
    
    % --- Column 3: Magnitude Spectra ---
    subplot(3, 3, (i-1)*3 + 3);
    plot(f_axis, P1_ref, 'k--', 'LineWidth', 1.2); hold on;
    plot(f_axis, P1_recon, 'b-', 'LineWidth', 1.1); hold off;
    title(sprintf('Magnitude Spectrum (f_s = %d Hz)', fs));
    xlabel('Frequency (Hz)'); ylabel('|X(f)|');
    xlim([0, 70]);
    grid on;
    if i == 1
        legend('Reference', 'Reconstructed', 'Location', 'northeast');
    end
end
