n = 2048;
a = 0.1;
NFFT = 2^nextpow2(n);
fs = 360;
test_sig = a * randn(1, n);
TEST = fftshift(fft(test_sig, NFFT))/n;
f = linspace(0, fs/2, NFFT/2 + 1);

subplot(4,1,1);
plot(f, abs(TEST(NFFT/2 : end)));

[B, A] = butter(4, 0.2);
y = filter(B, A, test_sig);
Y = fftshift(fft(y, NFFT))/n;
subplot(4,1,2);
plot(f, abs(Y(NFFT/2 : end)));

decimate = y(1 : 4 : end);
nfft_dec = 2^nextpow2(length(decimate));
DEC = fftshift(fft(decimate, nfft_dec))/length(decimate);
f1 = linspace(0, fs/8, nfft_dec/2 + 1);
subplot(4,1,3);
plot(f1, abs(DEC(nfft_dec/2 : end)));

reconstructed = zeros(1, n);
for i = 1 : length(decimate)
    reconstructed(4 * i) = decimate(i);
end
reconstructed = 4 * filter(B, A, reconstructed);

reconstruced_fft = fftshift(fft(reconstructed, NFFT))/n;
subplot(4,1,4);
plot(f, abs(reconstructed(NFFT/2 : end)));

figure;

subplot(2,1,1);
plot(1 : n , y);
xlim([1, n]);

subplot(2,1,2);
plot(1 : n, reconstructed);
xlim([1, n]);