tic
clc;
clear;
fsim = 1e8; % 
sim_time_step = 1/fsim ;
T = 1e-3;   %% Simulation time = 
freq_PNS = 1e6; % Nyquist Sampling frequency
Ts = 1/freq_PNS;

N = T*freq_PNS;
M = N/8 ; 
sparsity = 2;
%%        INPUT SIGNAL    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%f1 = 33000 ; f2 = 37000 ; f3 = 21000 ;
%f4 = 10000 ; f5 = 17000 ; f6 = 23000; f7 = 26000; f8 = 39000; f9 = 41000; f10 = 45000;
% f1 = 3300000000 ; f11 = 3700000000 ; f3 = 100000000 ;
% f4 = 1000000000 ; f5 = 1700000000 ; f6 = 2300000000; f7 = 2600000000; 
% f8 = 3900000000; f9 = 4100000000; f10 = 4500000000;
% 
t = 0:1/fsim:T;  
% x1 = sin(2*pi*f1*t);
% x2 = sin(2*pi*f11*t);
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
% x_sig = 6*x1 + x2 + 3*x3 + 8*x4 + 5*x5 + x6 + 8*x7 + x8 + 10*x9 + x10 ;
%x_sig = 5*x4 + 3*x2 ;
x_sig = generate_input(N,sparsity,freq_PNS,t)+offset;

noise = 2*randn(1,length(x_sig));
x = x_sig + noise;
disp('Input SNR');disp(snr(x_sig,noise));

ideal_x = x(1:fsim/freq_PNS:end-1);
ideal_x_sig = x_sig(1:fsim/freq_PNS:end-1);
NFFT = length(ideal_x);
ideal_X = fftshift(fft(ideal_x,NFFT)/NFFT);
ideal_X_sig = fftshift(fft(ideal_x_sig,NFFT)/NFFT);
%f = linspace(-freq_PNS/2,freq_PNS/2-freq_PNS/NFFT,NFFT);
f = linspace(0,freq_PNS/2 - freq_PNS/NFFT,NFFT/2);

subplot(3,1,1);
test = abs(ideal_X(NFFT/2 + 1 : NFFT));
plot(f,2*abs(ideal_X(NFFT/2+1 : NFFT)));
title(['Input signal FFT, N = ',num2str(N),', K = ',num2str(sparsity)]);
xlabel('f'); ylabel('|X(f)|');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RANDOM SEQUENCE GENERATION
% n = 7;  %% LFSR bit length
% ps = ones(1,n);
% ns = zeros(1,n); 
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
%%   RANDOM MODULATION  
prod = zeros(1,length(ideal_x));  %% Discrete Modulated Signal
i = 1;
while(i <= length(ideal_x))
    prod(i) = ideal_x(i)*rnd(1+mod(i,length(rnd)));
    i = i+1;
end
% PROD = fftshift(fft(prod)/length(prod));
% subplot(2,2,2);
% plot(f,abs(PROD(NFFT/2+1:NFFT)));
% title('Modulated signal spectrum');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTEGRATOR/LPF 
obs_samples = zeros(1,M);CF = N/M;
for i = 1:M
   obs_samples(i) = sum(prod((i-1)*CF+1 : CF*i)); 
end

% f2 = linspace(-freq_PNS/2,freq_PNS/2 - freq_PNS/M,M);
% OBS_SAMPLES = fftshift(fft(obs_samples)/M);
% subplot(2,2,4);
% plot(f2,abs(OBS_SAMPLES));
% title('Input to ADC (-M/N*fs/2 , M/N*fs/2)');xlabel('f'); 

obs_samples3 = prod(1:N/M:end-1);

Resolution = 6;
obs_samples1 = (max(abs(obs_samples))/(2^(Resolution-1)))*round((2^(Resolution-1))*obs_samples/max(abs(obs_samples))); %% Quantizing Integrator/Summation outputs
%obs_samples1 = round((2^(Resolution-1))*obs_samples/max(abs(obs_samples)));
SQNR = snr(obs_samples,obs_samples1 - obs_samples)
%% SIGNAL RECOVERY
theta = matrix_construct5(fsim,T,freq_PNS,rnd,N,M,Resolution);
%theta(:,400:600) = 0;
[sparse_indx,coeff] = dsp_3(theta,obs_samples1',sparsity);

% omp_out = algo_omp(sparsity, theta, obs_samples1');
% omp_out_id = (N - find(omp_out) + 1);
% amps = abs(nonzeros(omp_out));
% y = zeros(1, N/2);
% y(omp_out_id) = amps;
% subplot(3,1,2)
% plot(f, y);
% title(['Reconstructed signal FFT using OMP function']);
% xlabel('f'); ylabel('|X_{reconstructed}(f)|');

sparse_indx = (sparse_indx-1)*freq_PNS/N ;
%coeff = coeff(1:N/2) - (coeff(N:-1:N/2+1)); %%% !!!!!!!!!!!!!!!!
for i = 2:(N/2)
   coeff(i) = coeff(i) - coeff(N+2-i) ;
end

%disp(sparse_indx);
Sig_spectrum = 2*abs(ideal_X_sig(N/2+1:N));
Reconstructed_spectrum = abs(coeff(1:N/2));
Noise_spectrum = Reconstructed_spectrum-Sig_spectrum;
Output_SNR = 20*log10(norm(Sig_spectrum)/norm(Noise_spectrum))
%Output_SNR = 20*log10(norm(Reconstructed_spectrum)/norm(Noise_spectrum))

%f = linspace(-freq_PNS/2,freq_PNS/2-freq_PNS/NFFT,NFFT);
subplot(3,1,3);
plot(f,Reconstructed_spectrum);
title(['Reconstructed signal FFT using dsp_3 function (fs = ',num2str(freq_PNS),'  M = N/',num2str(N/M),') ']);
xlabel('f'); ylabel('|X''(f)|');

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
