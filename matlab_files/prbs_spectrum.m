%% 7, 10, 14 bit PRBS spectrum plot
clear;
% n = 7;
% ps = ones(1,n);
% ns = zeros(1,n);    
% %A = null(1,n);
% rnd = zeros(1, 2^n - 1);
% for i = 1:2^n-1
%     rnd(i) = ps(n)-xor(ps(n),1);
%     ns(2:n) = ps(1:n-1) ;
%     ns(1) = xor(ps(n),ps(1));
%     ps = ns ;
% end
range = [-1, 1];
band = [0 1];
prbs7 = idinput(2^7 - 1, 'prbs', band, range);
% prbs14 = (idinput(2^14 - 1, 'prbs', band, range))';
prbs11 = idinput(2^11 - 1, 'prbs', band, range);

n = 11;
ps = zeros(1,11);
ns = zeros(1,11); 
ckt_prbs = zeros(size(ps));
for i = 1 : 2^n - 1
    ckt_prbs(i) = ps(1);
    ns(2 : n) = ps(1 : n - 1);
    ns(1) = ~xor(ps(1), ps(3));
    ps = ns;
end
ckt_prbs = 2*ckt_prbs - 1;
prbs14 = ckt_prbs;
isequal(ckt_prbs', prbs11)
return;

%% 7b
% len = length(prbs7) - 1;
len = 2^10;
PRBS7 = fftshift(fft(prbs7,len));
f = linspace(-pi, pi - 2*pi/len,len);
subplot(3,1,1);
% plot(f/pi,abs(PRBS7).*abs(PRBS7));
plot(f/pi,abs(PRBS7).*abs(PRBS7), '-red', 'LineWidth', 1);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('|X_{fft}|^2', 'FontSize', 16);
xlabel('\times \pi', 'FontSize', 16);
grid on;
title('7 bit prbs spectrum', "FontSize", 18, "FontWeight", "bold");

% title('7 bit prbs spectrum');
% % ylim([0 35]);
% xlabel('\times \pi');
% ylabel('|X_{fft}|');

%% 11b
% len = length(prbs11) - 1;
% len = 2^10;
PRBS11 = fftshift(fft(prbs11,len));
f = linspace(-pi, pi - 2*pi/len,len);
subplot(3,1,2);
% plot(f/pi,abs(PRBS11).*abs(PRBS11));
plot(f/pi,abs(PRBS11).*abs(PRBS11), '-cyan', 'LineWidth', 1);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('|X_{fft}|^2', 'FontSize', 16);
xlabel('\times \pi', 'FontSize', 16);
grid on;
title('11 bit prbs spectrum', "FontSize", 18, "FontWeight", "bold");
%% 14b
% len = length(prbs14) - 1;
% len = 2^8; 
PRBS14 = fftshift(fft(prbs14,len));
f = linspace(-pi, pi - 2*pi/len,len);
subplot(3,1,3);
plot(f/pi,abs(PRBS14).*abs(PRBS14), '-black', 'LineWidth', 1);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('|X_{fft}|^2', 'FontSize', 16);
xlabel('\times \pi', 'FontSize', 16);
grid on;
title('14 bit prbs spectrum', "FontSize", 18, "FontWeight", "bold");
%% Periodograms
%{
figure;
[pxx7, w1] = periodogram(prbs7, [], length(prbs7) + 1);
subplot(3,1,1);
plot(w1, pxx7);

[pxx10, w2] = periodogram(prbs10, [], length(prbs10) + 1);
subplot(3,1,2);
plot(w2, pxx10);

[pxx14, w3] = periodogram(prbs14, [], length(prbs14) + 1);
subplot(3,1,3);
plot(w3, pxx14);
%}