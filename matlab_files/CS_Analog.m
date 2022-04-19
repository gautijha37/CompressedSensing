tic
clear;
fsim = 1e7; % 
sim_time_step = 1/fsim ;
T = 1e-3;   %% Simulation time frame
freq_PNS = 1e6; % Nyquist Sampling frequency
Ts = 1/freq_PNS;

N = T*freq_PNS;
M = N/8 ; 
sparsity = 20;
%% INPUT SIGNAL    
%f1 = 33000 ; f2 = 37000 ; f3 = 21000 ;
%f4 = 10000 ; f5 = 17000 ; f6 = 23000; f7 = 26000; f8 = 39000; f9 = 41000; f10 = 45000;
% f1 = 3300000000 ; f2 = 3700000000 ; f3 = 100000000 ;
% f4 = 1000000000 ; f5 = 1700000000 ; f6 = 2300000000; f7 = 2600000000; 
% f8 = 3900000000; f9 = 4100000000; f10 = 4500000000;
% 
t = 0:1/fsim:T;  
% x1 = sin(2*pi*f1*t);
% x2 = sin(2*pi*f2*t);
% x3 = sin(2*pi*f3*t);
% x4 = sin(2*pi*f4*t);
% x5 = sin(2*pi*f5*t);
% x6 = sin(2*pi*f6*t);
% x7 = sin(2*pi*f7*t);
% x8 = sin(2*pi*f8*t);
% x9 = sin(2*pi*f9*t);
% x10 = sin(2*pi*f10*t);
offset = 0*ones(1,length(t));
%x_sig = x1 ;
%x_sig = 6*x1 + x2 + 3*x3 + 8*x4 + 5*x5 + x6 + 8*x7 + x8 + 10*x9 + x10 ;
%x_sig = 5*x4 + 3*x2;
x_sig = generate_input(N,sparsity,freq_PNS,t)+offset;

noise = 2*randn(1,length(x_sig));
x = x_sig + noise;
Input_SNR = snr(x_sig,noise)

ideal_x = x(1:fsim/freq_PNS:end-1);
ideal_x_sig = x_sig(1:fsim/freq_PNS:end-1);
NFFT = length(ideal_x);
ideal_X = fftshift(fft(ideal_x,NFFT)/NFFT);
ideal_X_sig = fftshift(fft(ideal_x_sig,NFFT)/NFFT);
%f = linspace(-freq_PNS/2,freq_PNS/2-freq_PNS/NFFT,NFFT);
f = linspace(0,freq_PNS/2 - freq_PNS/NFFT,NFFT/2);

subplot(2,1,1);
plot(f,10*log10(2*abs(ideal_X(NFFT/2+1 : NFFT))));
title(['Input signal FFT, N = ',num2str(N),', K = ',num2str(sparsity)]);
xlabel('f'); ylabel('|X(f)| (dB)');
%% RANDOM SEQUENCE GENERATION 
% n = 7;  %% LFSR bit length
% ps = ones(1,n);
% ns = zeros(1,n);    
% %A = null(1,n);
% rnd = zeros(1,2^n - 1);
% for i = 1:2^n-1
%     rnd(i) = ps(n)-xor(ps(n),1);
%     ns(2:n) = ps(1:n-1) ;
%     ns(1) = xor(ps(n),ps(1));
%     ps = ns;
% end

range = [-1, 1];
band = [0 1];
prbs7 = idinput(2^7 - 1, 'prbs', band, range);
prbs10 = idinput(2^10 - 1, 'prbs', band, range);
rnd = idinput(2^14 - 1, 'prbs', band, range);

if(length(rnd) < T*freq_PNS)
    fprintf('PRBS length insufficient');
    return;
end

%%   RANDOM MODULATION 
prod = zeros(1,length(t));  %% Analog modulated Signal
i = 1; j = 1;
while(j <= length(rnd) && i <= length(t)) % length(rnd) must > T*freq_PNS
    for k = 1 : fsim/freq_PNS
        if(i > length(t))
            break;
        end
        prod(i) = x(i)*rnd(j);
        i = i+1;
    end
    j = j+1;
end

fprintf('PRBS length consumed = %d, i = %d', j, i);
%% INTEGRATOR/LPF 
% Integrates M consecutive windows independently in the simulation time
% frame
obs_samples = zeros(1,M);
i = 1; j = 1;
while (i < length(t)) % length(t) = T*fsim + 1
    obs_samples(j) = sum(prod(i : i + T*fsim/M - 1));
    j = j + 1;
    i = i + T*fsim/M ;
end
obs_samples = obs_samples*sim_time_step;
%% LOW RATE ADC 
Resolution = 6;
obs_samples1 = (max(abs(obs_samples))/(2^(Resolution-1)))*round((2^(Resolution-1))*obs_samples/max(abs(obs_samples))); %% Quantizing Integrator/Summation outputs
SQNR = snr(obs_samples,obs_samples - obs_samples1)
%% SIGNAL RECOVERY
theta = matrix_construct_analog(fsim,T,freq_PNS,rnd,N,M,Resolution);
[sparse_indx,coeff] = dsp_3(theta,obs_samples1',sparsity);
sparse_indx = (sparse_indx-1)*freq_PNS/N ;
%disp(sparse_indx);
Sig_spectrum = 2*abs(ideal_X_sig(N/2+1:N));
Reconstructed_spectrum = coeff(1:N/2);
Noise_spectrum = Reconstructed_spectrum-Sig_spectrum;
Output_SNR = 20*log10(norm(Sig_spectrum)/norm(Noise_spectrum))

subplot(2,1,2);
plot(f,10*log10(0.0001 + coeff(1:N/2)));
title(['Reconstructed signal FFT (fs = ',num2str(freq_PNS),'  M = N/',num2str(N/M),') ']);
xlabel('f'); ylabel('|X''(f)| (dB)');

%{
[XX,sparse_indx2] = OMP(10,obs_samples,phi);
sparse_indx2 = sparse_indx2*freq_PNS/N ;
disp(sparse_indx2);

XX = XX(1:N/2);
subplot(3,1,3);
plot(f,XX);
title(['Reconstructed signal FFT (fs = ',num2str(freq_PNS),'  M = N/',num2str(N/M),') ']);
xlabel('f'); ylabel('|X''(f)|');
%}
toc
