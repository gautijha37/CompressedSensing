%% Loading ECG Data
%{
    ECG data from : https://data.mendeley.com/datasets/7dybx7wyfn/3
    For research purposes, the ECG signals were obtained from the PhysioNet 
    service (http://www.physionet.org) from the MIT-BIH Arrhythmia database. 
    The created database with ECG signals is described below. 1) The ECG signals
    were from 45 patients: 19 female (age: 23-89) and 26 male (age: 32-89). 2) 
    The ECG signals contained 17 classes: normal sinus rhythm, pacemaker rhythm, 
    and 15 types of cardiac dysfunctions (for each of which at least 10 signal 
    fragments were collected). 3) All ECG signals were recorded at a sampling 
    frequency of 360 [Hz] and a gain of 200 [adu / mV]. 4) For the analysis, 
    1000, 10-second (3600 samples) fragments of the ECG signal (not overlapping) 
    were randomly selected. 5) Only signals derived from one lead, the MLII, were used. 6) Data are in mat format (Matlab).
%}
path = "C:\Users\gauti\OneDrive\Documents\MATLAB\MLII\1 NSR\100m (8).mat";
ecg_data = load(path);
% max(ecg_data.val)
% min(ecg_data.val)
sig = ecg_data.val/200; % in mV
sig = sig - mean(sig);
sig = sig(1 : 2048);
N = length(sig);
fs = 360; % Hz
% sparsity = find_ecg_peaks(sig); % Observed manually !!!
sparsity = 120;
t = linspace(0, (N - 1)/fs , N); %Time duration 10 sec.
subplot(2,2,1);
% plot(t, sig);
plot(t, sig, '-r', 'LineWidth', 1);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('Signal value (in mV)', 'FontSize', 16);
xlabel('Time (in sec)', 'FontSize', 16);
grid on;
title("Original ECG signal", "FontSize", 18, "FontWeight", "bold");

% xlabel('Time (in sec)');
% ylabel('Signal value (in mV)');
% title(['Original ECG Signal']);
%% Thresholding
W = 11; % Averaging window length
amp = max(sig) - min(sig);
th = 0.10 * amp; % Thresholding parameter
curr_avg = 0; % 0 for other averaging window
prev_num = 0;
threshold_sig = zeros(size(sig));
% for i = 1 : length(sig)
%     if((sig(i)) >= th + curr_avg) %% Ignored negative peaks, only considering R-R
%         threshold_sig(i) = sig(i);
%     else
%         threshold_sig(i) = 0;
%     end
%     
%     curr_avg = (curr_avg * W - prev_num + sig(i))/W;
%     if(i - W >= 0)
%         prev_num = sig(i - W + 1);
%     end
% end

curr_sum = sum(sig(1 : 1 + (W - 1)/2));
left_sum = 0;
i = 1;
rest_sig = zeros(size(sig));
mark = zeros(size(sig));
W2 = 11;
while( i <= length(sig))
    if((sig(i)) >= th + curr_sum/(2 * W) ) %% Ignored negative peaks, only considering R-R
        threshold_sig(i - (W - 1)/2 : i + (W - 1)/2) = sig(i - (W - 1)/2 : i + (W - 1)/2);
%         rest_sig(i) = left_sum/W2;
        mark(i - (W - 1)/2 : i + (W - 1)/2) = ones(1, W);
    end
    
    if(i + (W - 1)/2 + 1 <= length(sig))
        curr_sum = curr_sum + sig(i + (W - 1)/2 + 1);
    end
    
    if(i - (W - 1)/2 >= 1)
        curr_sum = curr_sum - sig(i - (W - 1)/2);
    end
    i = i + 1;
end

for i = 1 : length(sig)
    if(mark(i) == 0)
        rest_sig(i) = sig(i);
    else
        rest_sig(i) = left_sum/W2;
    end
    
    if(i >= W2 + 1)
        left_sum = left_sum - rest_sig(i - W2);
    end
    
    left_sum = left_sum + rest_sig(i);
end

%% Adding AWGN
SNR = 30; %dB
noise_rms = rms(threshold_sig)/(10^(SNR/20));
noise = noise_rms * randn(size(threshold_sig));
% rms(noise)
noisy_sig = noise + threshold_sig;

subplot(2,2,2);
% plot(t, noisy_sig);
plot(t, threshold_sig, '-cyan', 'LineWidth', 1);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('Signal value (in mV)', 'FontSize', 16);
xlabel('Time (in sec)', 'FontSize', 16);
grid on;
title(['Thresholded ECG signal with W = ', num2str(W), ', th = ', num2str(100*th/amp), ' %'], "FontSize", 18, "FontWeight", "bold");

% title(['Thresholded noisy ECG signal with W = ', num2str(W), ', th = ', num2str(100*th/amp), ' %']);
% xlabel('Time (in sec)');
% ylabel('Signal value (in mV)');
%% Compression
CF = 4;
M = N/CF;

n = 11; % PRBS bit length
range = [-1, 1];
band = [0 1];
prbs11 = idinput(2^11 - 1, 'prbs', band, range);
% PRBS_mat = construct_PRBS_mat(N, M, 64, prbs11');

PRBS_mat = construct_full_PRBS_mat(M, N, n);
gaussian_mat = round(2 * (randn(M, N) > 0) - 1, 5); % digits of precision in gaussian matrix
theta = PRBS_mat; % Psi = eye(N, N)
% theta = gaussian_mat;

observed_sig = theta * threshold_sig';
% subplot(2,2,3);
% % plot(1 : M, observed_sig);
% plot(1 : M, observed_sig, '-black', 'LineWidth', 1);
% set(gca, 'FontSize', 8, 'FontWeight', 'bold');
% ylabel('Compressed Signal', 'FontSize', 16);
% xlabel('Sample No.', 'FontSize', 16);
% xlim([1 512]);
% grid on;
% title(['Compressed ECG signal with CF = ', num2str(CF)], "FontSize", 18, "FontWeight", "bold");

% title(['Compressed ECG signal with CF = ', num2str(CF)]);
% xlabel('Sample No.');
% ylabel('Compressed Signal');

%% Frequency sparse part of signal
subplot(2,2,3);
% plot(1 : M, observed_sig);
plot(t , rest_sig, '-black', 'LineWidth', 1);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('Residual Signal (in mV)', 'FontSize', 16);
xlabel('time (in sec)', 'FontSize', 16);
% xlim([1 512]);
grid on;
title(['Residual ECG signal with CF = ', num2str(CF)], "FontSize", 18, "FontWeight", "bold");
%% Reconstruction
% Testing for Synthetic time domain sparse signal
% rand_id = randperm(length(sig), sparsity);
% rand_vals = 0.1 * max(sig) * rand(1, sparsity) + 0.9 * max(sig);
% test_sig = zeros(size(sig));
% test_sig(rand_id) = rand_vals;
% test_sig = test_sig + noise;
% observed_test_sig = theta * test_sig';

% subplot(3,1,1);
% plot(t, test_sig);
% xlabel('Time in sec');
% ylabel('Signal amplitude in mV');
% title('TEST ECG signal');

% subplot(3,1,2);
% plot(1 : M, observed_test_sig);
% xlabel('Observed Sample No.');
% ylabel('Measured test signal in mV');
% title('MEASURED TEST ECG signal');


% [sparse_indx,coeff] = dsp_3(theta,observed_sig,sparsity);
% sparse_indx
% find(test_sig)
% cc = coeff;

% sparse_indx = (sparse_indx-1)*freq_PNS/N ;
% for i = 2:(N/2)
%    coeff(i) = coeff(i) - coeff(N+2-i) ;
% end
% Reconstructed_spectrum = abs(coeff(1:N/2));
% scaling_factor = 40; % Seen from observation of FFT of Input sine of amplitude Vamp
subplot(2,2,4);
y =  (algo_omp(sparsity, theta, observed_sig))';
% for i = 1 : length(y)
%     if(abs(y(i)) < amp/32)
%         y(i) = 0;
%     end
% end
RSNR = snr(threshold_sig, threshold_sig - y)
% y = run_omp(t, theta, observed_sig);
% plot(t, y);
plot(t, y, '-magenta', 'LineWidth', 1);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('Signal amplitude in mV', 'FontSize', 16);
xlabel('Time in sec', 'FontSize', 16);
grid on;
title('Reconstructed ECG signal', "FontSize", 18, "FontWeight", "bold");

% xlabel('Time in sec');
% ylabel('Signal amplitude in mV');
% title('Reconstructed ECG signal');

%{
figure;
subplot(2, 1, 1);
plot(t, sig, '-r', 'LineWidth', 1);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('Signal value (in mV)', 'FontSize', 16);
xlabel('Time (in sec)', 'FontSize', 16);
grid on;
title("Original ECG signal", "FontSize", 18, "FontWeight", "bold");

subplot(2,1,2);
plot(t, threshold_sig, '-b', 'LineWidth', 1);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('Signal value (in mV)', 'FontSize', 16);
xlabel('Time (in sec)', 'FontSize', 16);
grid on;
title(['Thresholded ECG signal with W = ', num2str(W), ', th = ', num2str(100*th/amp), ' %'], "FontSize", 18, "FontWeight", "bold");
%}
% title(['Thresholded noisy ECG signal with W = ', num2str(W), ', th = ', num2str(100*th/amp), ' %']);
% xlabel('Time (in sec)');
% ylabel('Signal value (in mV)');
% 
% plot(t, BER_h_actual, '-*g', 'MarkerSize', 16, 'LineWidth', 2);
% plot(SNR_range, BER_no_h, '-*b', 'MarkerSize', 16, 'LineWidth', 2);
% set(gca, 'FontSize', 16, 'FontWeight', 'bold');
% ylabel('BER', 'FontSize', 20), xlabel('Eb/N0 (dB)', 'FontSize', 20);
% grid on;
% legend({'Equalized (Using estimated H) ', 'Equalized (Using actual H)', 'No Equalization'}, 'FontSize', 14, 'Location', 'best');
% title("BER curves for Rayleigh Channel", "FontSize", 24, "FontWeight", "bold"); hold off;


%% Processing residue signal
figure;
subplot(2,1,1);
% plot(1 : M, observed_sig);
plot(t , rest_sig, '-black', 'LineWidth', 1);
set(gca, 'FontSize', 5, 'FontWeight', 'bold');
ylabel('Residual Signal (in mV)', 'FontSize', 10);
xlabel('time (in sec)', 'FontSize', 10);
% xlim([1 512]);
grid on;
title(['Residual ECG signal with CF = ', num2str(CF), ' using modified thresholding algorithm'], "FontSize", 18, "FontWeight", "bold");

NFFT = length(rest_sig);
REST_SIG = fftshift(fft(rest_sig,NFFT)/NFFT);
% ideal_X_sig = fftshift(fft(ideal_x_sig,NFFT)/NFFT);
%f = linspace(-freq_PNS/2,freq_PNS/2-freq_PNS/NFFT,NFFT);
f = linspace(-fs/2,fs/2 - fs/NFFT,NFFT);

% subplot(3,1,2);
% plot(t , sig - threshold_sig, '-black', 'LineWidth', 1);
% set(gca, 'FontSize', 5, 'FontWeight', 'bold');
% ylabel('Original Residue', 'FontSize', 10);
% xlabel('time in sec', 'FontSize', 10);
% % xlim([1 512]);
% grid on;
% title(['Residual ECG signal'], "FontSize", 10, "FontWeight", "bold");

[B, A] = butter(4, 0.2);
filtered_residue = filtfilt(B, A, rest_sig);
downsample = filtered_residue(1 : 4 : end);

reconstructed = zeros(1, N);
for i = 1 : length(downsample)
    reconstructed(4 * i) = downsample(i);
end
reconstructed = 4 * filtfilt(B, A, reconstructed);

subplot(2,1,2);
% f = linspace(0,fs/2 - fs/N, N/2);
plot(t,reconstructed, '-black', 'LineWidth', 1);
set(gca, 'FontSize', 5, 'FontWeight', 'bold');
title(['Received residue signal after upsampling'], "FontSize", 18, "FontWeight", "bold");
ylabel('Reconstructed residue', 'FontSize', 10);
xlabel('Time (in sec)', 'FontSize', 10);

figure;
subplot(2,1,1);
% plot(t, sig);
plot(t, sig, '-black', 'LineWidth', 1);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('Signal value (in mV)', 'FontSize', 16);
xlabel('Time (in sec)', 'FontSize', 16);
grid on;
title("Original ECG signal", "FontSize", 18, "FontWeight", "bold");

subplot(2,1,2);
% plot(t, sig);
final_sig = y;
for i = 1 : length(sig)
    if(final_sig(i) == 0)
        final_sig(i) = reconstructed(i);
    end
end
plot(t, final_sig, '-black', 'LineWidth', 1);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('Signal value (in mV)', 'FontSize', 16);
xlabel('Time (in sec)', 'FontSize', 16);
grid on;
title("Final reconstructed signal", "FontSize", 18, "FontWeight", "bold");

final_snr = snr(sig, final_sig - sig)
% % Project rest_sig on smaller dimension M
% range = [-1, 1];
% band = [0 1];
% rnd = idinput(2^14 - 1, 'prbs', band, range);
% 
% prod = zeros(1,length(rest_sig));  %% Discrete Modulated Signal
% i = 1;
% while(i <= length(rest_sig))
%     prod(i) = rest_sig(i)*rnd(1+mod(i,length(rnd)));
%     i = i+1;
% end
% 
% obs_samples = zeros(1,M);
% for i = 1:M
%    obs_samples(i) = sum(prod((i-1)*CF+1 : CF*i)); 
% end
% 
% Resolution = 6;
% % T = t(end);
% % fsim = 1e8;
% theta_freq = matrix_construct6(rest_sig, M, rnd, fs);
% sparsity = 4;
% [sparse_indx,coeff] = dsp_3(theta,obs_samples',sparsity);
% 
% sparse_indx = (sparse_indx-1)*fs/N ;
% %coeff = coeff(1:N/2) - (coeff(N:-1:N/2+1)); %%% !!!!!!!!!!!!!!!!
% for i = 2:(N/2)
%    coeff(i) = coeff(i) - coeff(N+2-i) ;
% end
% 
% Reconstructed_spectrum = abs(coeff(1:N/2));
% 
% subplot(3,1,3);
% f = linspace(0,fs/2 - fs/N, N/2);
% plot(f,Reconstructed_spectrum, '-black', 'LineWidth', 1);
% set(gca, 'FontSize', 5, 'FontWeight', 'bold');
% title(['Residual signal FFT using dsp_3 function (fs = ',num2str(fs)], "FontSize", 10, "FontWeight", "bold");
% ylabel('|X''(f)|', 'FontSize', 10);
% xlabel('Freq (Hz)', 'FontSize', 10);


