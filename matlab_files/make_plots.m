s1 = load('plot.mat');
s2 = load('plot2.mat');


subplot(2,2,1);
plot(s1.f, abs(s1.VIN), '-red', 'LineWidth', 3);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('|X_{in}(f)|', 'FontSize', 16);
xlabel('f', 'FontSize', 16);
xlim([0 12.5e6]);
grid on;
title('Input Spectrum', "FontSize", 18, "FontWeight", "bold");

subplot(2,2,2);
plot(s2.f, abs(s2.VIN), '-red', 'LineWidth', 3);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('|X_{in}(f)|', 'FontSize', 16);
xlabel('f', 'FontSize', 16);
xlim([0 12.5e6]);
grid on;
title('Input Spectrum', "FontSize", 18, "FontWeight", "bold");

subplot(2,2,3);
plot(s1.f, s1.scaling_factor*s1.Reconstructed_spectrum, '-black', 'LineWidth', 3);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('|X(f)|', 'FontSize', 16);
xlabel('f', 'FontSize', 16);
xlim([0 12.5e6]);
grid on;
title('Reconstructed Spectrum', "FontSize", 18, "FontWeight", "bold");

subplot(2,2,4);
plot(s2.f, s2.scaling_factor*s2.Reconstructed_spectrum, '-black', 'LineWidth', 3);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('|X(f)|', 'FontSize', 16);
xlabel('f', 'FontSize', 16);
xlim([0 12.5e6]);
grid on;
title('Reconstructed Spectrum', "FontSize", 18, "FontWeight", "bold");
