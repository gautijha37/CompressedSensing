clc
clear all
%% netlist path %%
%"C:\Users\gauti\Documents\MATLAB\cmi_7b_2tone.net"
%% .exe path %%
%"C:\Program Files\LTC\LTspiceXVII\scad3.exe"

%% My netlist
netlist = ['C:\Users\gauti\OneDrive\Documents\MATLAB\cmi_7b_2tone.net'];

input_amp1 = 0.4; % Input sine amp in volt
input_freq1 = 0.46; % fin = input_freq * fs

input_amp2 = 0.3; % Input sine amp in volt
input_freq2 = 0.31; % fin = input_freq * fs

data = fileread('netlist_7b_2tone.txt');
txt_data2 = strsplit(data,'\n');
str = string(txt_data2);
string_to_write = ['.param X1 = ', num2str(input_freq1), ' X2 = ', num2str(input_freq2), newline];
string_for_amp = ['.param vamp1 = ', num2str(input_amp1), ' vamp2 = ', num2str(input_amp2), newline];
str(end : end + 2) = str(end - 2:end);
str(end-3) = string_to_write;
str(end-4) = string_for_amp;

fid = fopen(netlist,'w+');
fprintf(fid, '%s\n', str);
fid = fclose(fid);

fileName = 'cmi_7b_2tone.raw';
dos('LTSpice_call_7bit_2tone.bat');
pause(10) ;
raw_data = LTspice2Matlab(fileName);
dos('LTSpice_end.bat');

%% Extracting signal samples %%
VDC = 0.9; T = 40e-9; Tp = T/2 ;
Vamp = 0.4; CF = 4;Tp1 = Tp/2;
N = round(raw_data.time_vect(end)/T);
M = N/CF; G = 1/(10e3); C = 5e-12;
freq_PNS = 1/T ; clk_rise_time = 500e-12;

id_vin = find(contains(raw_data.variable_name_list,'V(vin)'));
id_vp1 = find(contains(raw_data.variable_name_list,'V(vp1)'));
id_vp2 = find(contains(raw_data.variable_name_list,'V(vp2)'));
id_v_snh = find(contains(raw_data.variable_name_list,'V(v_snh)'));
id_prbs = find(contains(raw_data.variable_name_list,'V(prbs_out)'));

vp1 = raw_data.variable_mat(id_vp1, :);
vp2 = raw_data.variable_mat(id_vp2, :);
vin = raw_data.variable_mat(id_vin, :);
v_snh = raw_data.variable_mat(id_v_snh, :);
prbs_out = raw_data.variable_mat(id_prbs, :);
% plot(raw_data.time_vect,vin);

samples_vp = zeros(1, N); % Combined outputs of vp1 and vp2 after each clock cycle
i = T + Tp/2; j = 1;
while(i < raw_data.time_vect(end))
    for ii = 1 : CF
        if(i < raw_data.time_vect(end))
            sample_id = find(abs(raw_data.time_vect - i) == min(abs(raw_data.time_vect - i))) ;
            samples_vp(j) = vp1(sample_id);
            i = i + T; j = j + 1;
        end
    end
    for ii = 1 : CF
        if(i < raw_data.time_vect(end))
            sample_id = find(abs(raw_data.time_vect - i) == min(abs(raw_data.time_vect - i))) ;
            samples_vp(j) = vp2(sample_id);
            i = i + T; j = j + 1;
        end
    end
end
samples_vp_copy = samples_vp;
samples_vp = samples_vp - VDC;

ckt_prbs = zeros(1, N);
i = Tp; j = 1;
while(i < raw_data.time_vect(end))
    sample_id = find(abs(raw_data.time_vect - i) == min(abs(raw_data.time_vect - i))) ;
    ckt_prbs(j) = prbs_out(sample_id);
    i = i + T; j = j + 1;
end
ckt_prbs = round(ckt_prbs/1.8);
ckt_prbs = 2*ckt_prbs - 1;
    
obs_samples_vp1 = zeros(1, M); % Output after 4 clock cycles at vp1
i = CF*T + Tp/2; j = 1;
while(i < raw_data.time_vect(end))
   sample_id = find(abs(raw_data.time_vect - i) == min(abs(raw_data.time_vect - i))) ;
    obs_samples_vp1(j) = vp1(sample_id);
   i = i + 2*CF*T; j = j + 1;
end

obs_samples_vp2 = zeros(1, M); % Output after 4 clock cycles at Vp2
i = 2*CF*T + Tp/2; j = 1;
while(i < raw_data.time_vect(end))
   sample_id = find(abs(raw_data.time_vect - i) == min(abs(raw_data.time_vect - i))) ;
    obs_samples_vp2(j) = vp2(sample_id);
   i = i + 2*CF*T; j = j + 1;
end

obs_samples_vp = zeros(1, M); % Combined outputs of vp1 and vp2
i = 1; j = 1;
while (i <= M)
    obs_samples_vp(i) = obs_samples_vp1(j);
    obs_samples_vp(i + 1) = obs_samples_vp2(j);
    j = j + 1;
    i = i + 2;
end
obs_samples_vp_copy = obs_samples_vp;

obs_samples_vin = zeros(1, N); j = 1;
i = Tp;
while(i < raw_data.time_vect(end))
    sample_id = find(abs(raw_data.time_vect - i) == min(abs(raw_data.time_vect - i))) ;
    obs_samples_vin(j) = vin(sample_id);
    j = j + 1; i = i + T;
end
obs_samples_vin = obs_samples_vin - VDC;

obs_samples_v_snh = zeros(1, N); j = 1;
i = 1.5*Tp;
while(i < raw_data.time_vect(end))
    sample_id = find(abs(raw_data.time_vect - i) == min(abs(raw_data.time_vect - i))) ;
    obs_samples_v_snh(j) = v_snh(sample_id);
    j = j + 1; i = i + T;
end
t1 = obs_samples_v_snh;

f = linspace(0,freq_PNS/2 - freq_PNS/N, N/2);
subplot(2,1,1);
VIN = fftshift(fft(obs_samples_vin,N));
VIN = VIN(N/2 + 1 : N);
plot(f, abs(VIN), '-red', 'LineWidth', 3);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('|X_{in}(f)|', 'FontSize', 16);
xlabel('f', 'FontSize', 16);
grid on;
title('Input Spectrum', "FontSize", 18, "FontWeight", "bold");
% subplot(2,1,2);
% plot(f, 20*log10(0.001 + abs(fftshift(fft(obs_samples_v_snh,N)/N))));

%% Theta construction %%

sample_times_clk = Tp : T : raw_data.time_vect(end);
% [PRBS_mat, rnd_data] = construct_PRBS_mat(N, CF);
[PRBS_mat, rnd] = construct_PRBS_mat(N, CF);
theta = (Vamp*G*Tp/C)*construct_theta(PRBS_mat, freq_PNS,sample_times_clk, VDC);
% corr_mat_theta = correl(theta);
% struc_data = load('cmi53.mat','theta_data');
% theta_data = struc_data.theta_data;
% theta_data = construct_theta_data(N, CF, fileName, txt_data2, netlist, Tp, T, VDC, id_vp);
% theta_data = theta_data - repmat(PRBS_mat * (G*Tp/C*VDC*(ones(N, 1))), [1, N]);

rnd = [rnd, rnd, rnd, rnd]; % Making rnd length 4x.
ideal_vp_vals = zeros(1, N);
i = 1;
while(i <= N)
    for j = 1 : CF
        ideal_vp_vals(i + j - 1) = rnd(i : i + j - 1) * obs_samples_v_snh(i : i + j - 1)';
    end
    i = i + CF;
end
ideal_vp_vals = G*Tp/C * ideal_vp_vals;
error1 = samples_vp - ideal_vp_vals ;
% plot(1 : N, error1);

% theta_data_id = zeros(size(theta_data));
% norm_vals = zeros(1, N);
% for i = 1:N
%     norm_vals(i) = norm(theta_data(:, i));
%     theta_data_id(:, i) = theta_data(:, i)/norm_vals(i);
% end    
% corr_mat_theta_data = correl(theta_data);
% corr_mat_theta_data2 = correl(theta_data_id);

obs_samples_vp = obs_samples_vp - VDC - (PRBS_mat * (G*VDC*Tp/C*(ones(N, 1))))';
test_vin = PRBS_mat * obs_samples_vin';
test_v_snh = PRBS_mat * obs_samples_v_snh';
t1 = (G*Tp/C)*PRBS_mat * t1';
t1_final = t1 - PRBS_mat * (Tp/C*G*VDC*(ones(N, 1)));
%% Signal Reconstruction %%
% new_obs_samples_vp = obs_samples_vp + (PRBS_mat * (G*VDC*(ones(N, 1))))'
% figure;

[sparse_indx,coeff] = dsp_3(theta,obs_samples_vp',2);
% [sparse_indx,coeff] = dsp_4(theta_data_id, theta_data ,obs_samples_vp',1);
% [sparse_indx,coeff] = dsp_3(theta,(G*Tp/C)*test_vin,1);
% [sparse_indx,coeff] = dsp_3(theta,(G*Tp/C)*test_v_snh,1);
% [sparse_indx,coeff] = dsp_3(theta,t1_final,1);
sparse_indx
cc = coeff;
% aa = t1 - t2';
sparse_indx = (sparse_indx-1)*freq_PNS/N ;
for i = 2:(N/2)
   coeff(i) = coeff(i) - coeff(N+2-i) ;
end
Reconstructed_spectrum = abs(coeff(1:N/2));
scaling_factor = 40; % Seen from observation of FFT of Input sine of amplitude Vamp
subplot(2,1,2);
plot(f, scaling_factor*Reconstructed_spectrum, '-black', 'LineWidth', 3);
set(gca, 'FontSize', 8, 'FontWeight', 'bold');
ylabel('|X(f)|', 'FontSize', 16);
xlabel('f', 'FontSize', 16);
grid on;
title('Reconstructed Spectrum', "FontSize", 18, "FontWeight", "bold");
