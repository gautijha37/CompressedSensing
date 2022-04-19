clc
clear all
%% netlist path %%
%"C:\Users\Gautam Jha\Documents\LTspiceXVII\Current_Mode_Integrator_3.net"
%% .exe path %%
%"C:\Program Files\LTC\LTspiceXVII\scad3.exe"

%% My netlist
%{
* C:\Users\Gautam Jha\Documents\LTspiceXVII\Current_Mode_Integrator_3.asc
A1 N002 VSS VCLK VSS VSS ~PRBS_OUT PRBS_OUT VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A2 PRBS_OUT VSS VCLK VSS VSS VSS N003 VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A3 N003 VSS VCLK VSS VSS VSS N004 VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A4 N004 VSS VCLK VSS VSS VSS N005 VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A5 N005 VSS VCLK VSS VSS VSS N006 VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A6 N006 VSS VCLK VSS VSS VSS N007 VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A7 N007 VSS VCLK VSS VSS VSS N001 VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A8 N001 PRBS_OUT VSS VSS VSS N002 VSS VSS XOR Vhigh = 1.8 Vlow = 0
V1 VCLK 0 PULSE({VSS} {VDD} {Tp} 50p 50p {Tp} {T})
V2 VDD 0 {VDD}
V3 VSS 0 {VSS}
M1 N013 ~CHARGE Vp VDD CMOSP l=0.5u w=6u
M2 Vp DISCHARGE N015 VSS CMOSN l=.5u w=6u
C1 Vp 0 5p
XU1 Vpc Vp Vpc opamp Aol=500 GBW=300Meg
M3 N014 Vbn1 VSS VSS CMOSN l={Ln1} w={Wn1} m={fn1}
M4 VDD Vbp1 Vbp1 VDD CMOSP l={Lp1} w={Wp1} m={fp1}
M5 N016 Vbn1 VSS VSS CMOSN l={Ln1} w={Wn1} m={fn1}
M6 Vpc CHARGE N013 VDD CMOSP l=.5u w=6u
M7 N015 ~DISCHARGE Vpc VSS CMOSN l=.5u w=6u
M8 VDD Vbp1 N010 VDD CMOSP l={Lp1} w={Wp1} m={fp1}
M9 Vbp2 Vbn2 N014 VSS CMOSN l={Ln2} w={Wn2} m={fn2}
M10 Vbp1 Vbp2 Vbp2 VDD CMOSP l={Lp2} w={Wp2} m={fp2}
M11 N010 Vbp2 N013 VDD CMOSP l={Lp2} w={Wp2} m={fp2}
M12 N015 Vbn2 N016 VSS CMOSN l={Ln2} w={Wn2} m={fn2}
M13 VDD N011 N012 VDD CMOSP l=1u w=6u
R1 N012 VSS 10K
M14 VDD N011 Vbn2 VDD CMOSP l=1u w=6u
M15 Vbn1 Vbn1 VSS VSS CMOSN l={Ln1} w={Wn1} m={fn1}
XU2 V_SnH N012 N011 opamp Aol=500 GBW=200Meg
M16 Vbn2 Vbn2 Vbn1 VSS CMOSN l={Ln2} w={Wn2} m={fn2}
A9 VCLK3 PRBS_OUT VDD VDD VDD 0 CHARGE 0 AND Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A10 VCLK3 ~PRBS_OUT VDD VDD VDD 0 DISCHARGE 0 AND Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A11 CHARGE 0 0 0 0 ~CHARGE 0 0 BUF Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A12 DISCHARGE 0 0 0 0 ~DISCHARGE 0 0 BUF Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A13 VCLK 0 0 0 0 ~VCLK 0 0 BUF Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
V4 VCLK1 0 PULSE({VSS} {VDD} 0 50p 50p {0.425*T} {T})
A14 VCLK1 0 0 0 0 ~VCLK1 0 0 BUF Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
V5 VCLK2 0 PULSE({VSS} {VDD} 0 50p 50p {0.375*T} {T})
V6 VCLK3 0 PULSE({VSS} {VDD} {0.475*T} 50p 50p {0.475*T} {T})
V7 VDC 0 0.9
M17 Va VCLK1 Vin VSS CMOSN l=1u w=4u
M18 Vin ~VCLK1 Va VDD CMOSP l=1u w=4u
A15 VCLK2 0 0 0 0 ~VCLK2 0 0 BUF Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A16 VCLK3 0 0 0 0 ~VCLK3 0 0 BUF Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
C2 Vb Va 0.1p
XU3 Vb VDC V_SnH opamp Aol=500 GBW=300mega
M19 V_SnH VCLK2 Vb VSS CMOSN l=1u w=4u
M20 Vb ~VCLK2 V_SnH VDD CMOSP l=1u w=4u
M21 V_SnH VCLK3 Va VSS CMOSN l=1u w=4u
M22 Va ~VCLK3 V_SnH VDD CMOSP l=1u w=4u
V8 Vin VDC SINE(0 0.3 {0.125/T})
M23 VDC RESET Vp VSS CMOSN l=1u w=4u m=3
M24 Vp ~RESET VDC VDD CMOSP l=1u w=4u m=3
A17 N008 0 VCLK1 VSS VSS N008 Q0 0 DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A18 N009 0 Q0 VSS VSS N009 Q1 0 DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
A19 0 Q1 Q0 VCLK1 VDD ~RESET RESET 0 AND Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF
.model NMOS NMOS
.model PMOS PMOS
.lib C:\Users\Gautam Jha\Documents\LTspiceXVII\lib\cmp\standard.mos
.param VDD=1.8V VSS=0V td_DFF=1p T=30n Tp=15n tff_rise=50p
.ic V(Vp) = 0.9 V(Vpc) = 0.9
.include tsmc018.txt
.tran {8*T}
.LIB opamp.sub
.param Wn1=4u Ln1=0.25u Wn2=4u Ln2=.25u Wp1=20u Lp1=0.25u Wp2=20u Lp2=0.25u
.param fn1=8 fn2=8 fp1=8 fp2=8
.ic V(Va) = 0.9 V(Vb) = 0.9 V(V_SnH) = 0.9
.backanno
.end

%}

netlist = ['C:\Users\Gautam Jha\Documents\LTspiceXVII\Current_Mode_Integrator_3.net'];
%{
code = ['* C:\\Users\\Gautam Jha\\Documents\\LTspiceXVII\\MATLAB_Ltspice\\Current_Mode_Integrator_3.asc                            \r\n'...
'A1 N002 VSS VCLK VSS VSS ~PRBS_OUT PRBS_OUT VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A2 PRBS_OUT VSS VCLK VSS VSS VSS N003 VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A3 N003 VSS VCLK VSS VSS VSS N004 VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A4 N004 VSS VCLK VSS VSS VSS N005 VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A5 N005 VSS VCLK VSS VSS VSS N006 VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A6 N006 VSS VCLK VSS VSS VSS N007 VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A7 N007 VSS VCLK VSS VSS VSS N001 VSS DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A8 N001 PRBS_OUT VSS VSS VSS N002 VSS VSS XOR Vhigh = 1.8 Vlow = 0\r\n'...
'V1 VCLK 0 PULSE({VSS} {VDD} {Tp} 50p 50p {Tp} {T})\r\n'...
'V2 VDD 0 {VDD}\r\n'...
'V3 VSS 0 {VSS}\r\n'...
'M1 N013 ~CHARGE Vp VDD CMOSP l=0.5u w=6u\r\n'...
'M2 Vp DISCHARGE N015 VSS CMOSN l=.5u w=6u\r\n'...
'C1 Vp 0 5p\r\n'...
'XU1 Vpc Vp Vpc opamp Aol=500 GBW=300Meg\r\n'...
'M3 N014 Vbn1 VSS VSS CMOSN l={Ln1} w={Wn1} m={fn1}\r\n'...
'M4 VDD Vbp1 Vbp1 VDD CMOSP l={Lp1} w={Wp1} m={fp1}\r\n'...
'M5 N016 Vbn1 VSS VSS CMOSN l={Ln1} w={Wn1} m={fn1}\r\n'...
'M6 Vpc CHARGE N013 VDD CMOSP l=.5u w=6u\r\n'...
'M7 N015 ~DISCHARGE Vpc VSS CMOSN l=.5u w=6u\r\n'...
'M8 VDD Vbp1 N010 VDD CMOSP l={Lp1} w={Wp1} m={fp1}\r\n'...
'M9 Vbp2 Vbn2 N014 VSS CMOSN l={Ln2} w={Wn2} m={fn2}\r\n'...
'M10 Vbp1 Vbp2 Vbp2 VDD CMOSP l={Lp2} w={Wp2} m={fp2}\r\n'...
'M11 N010 Vbp2 N013 VDD CMOSP l={Lp2} w={Wp2} m={fp2}\r\n'...
'M12 N015 Vbn2 N016 VSS CMOSN l={Ln2} w={Wn2} m={fn2}\r\n'...
'M13 VDD N011 N012 VDD CMOSP l=1u w=6u\r\n'...
'R1 N012 VSS 10K\r\n'...
'M14 VDD N011 Vbn2 VDD CMOSP l=1u w=6u\r\n'...
'M15 Vbn1 Vbn1 VSS VSS CMOSN l={Ln1} w={Wn1} m={fn1}\r\n'...
'XU2 V_SnH N012 N011 opamp Aol=500 GBW=200Meg\r\n'...
'M16 Vbn2 Vbn2 Vbn1 VSS CMOSN l={Ln2} w={Wn2} m={fn2}\r\n'...
'A9 VCLK3 PRBS_OUT VDD VDD VDD 0 CHARGE 0 AND Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A10 VCLK3 ~PRBS_OUT VDD VDD VDD 0 DISCHARGE 0 AND Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A11 CHARGE 0 0 0 0 ~CHARGE 0 0 BUF Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A12 DISCHARGE 0 0 0 0 ~DISCHARGE 0 0 BUF Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A13 VCLK 0 0 0 0 ~VCLK 0 0 BUF Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'V4 VCLK1 0 PULSE({VSS} {VDD} 0 50p 50p {0.425*T} {T})\r\n'...
'A14 VCLK1 0 0 0 0 ~VCLK1 0 0 BUF Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'V5 VCLK2 0 PULSE({VSS} {VDD} 0 50p 50p {0.375*T} {T})\r\n'...
'V6 VCLK3 0 PULSE({VSS} {VDD} {0.475*T} 50p 50p {0.475*T} {T})\r\n'...
'V7 VDC 0 0.9\r\n'...
'M17 Va VCLK1 Vin VSS CMOSN l=1u w=4u\r\n'...
'M18 Vin ~VCLK1 Va VDD CMOSP l=1u w=4u\r\n'...
'A15 VCLK2 0 0 0 0 ~VCLK2 0 0 BUF Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A16 VCLK3 0 0 0 0 ~VCLK3 0 0 BUF Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'C2 Vb Va 0.1p\r\n'...
'XU3 Vb VDC V_SnH opamp Aol=500 GBW=300mega\r\n'...
'M19 V_SnH VCLK2 Vb VSS CMOSN l=1u w=4u\r\n'...
'M20 Vb ~VCLK2 V_SnH VDD CMOSP l=1u w=4u\r\n'...
'M21 V_SnH VCLK3 Va VSS CMOSN l=1u w=4u\r\n'...
'M22 Va ~VCLK3 V_SnH VDD CMOSP l=1u w=4u\r\n'...
'V8 Vin VDC SINE(0 0.3 {0.125/T})\r\n'...
'M23 VDC RESET Vp VSS CMOSN l=1u w=4u m=3\r\n'...
'M24 Vp ~RESET VDC VDD CMOSP l=1u w=4u m=3\r\n'...
'A17 N008 0 VCLK1 VSS VSS N008 Q0 0 DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A18 N009 0 Q0 VSS VSS N009 Q1 0 DFLOP Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'...
'A19 0 Q1 Q0 VCLK1 VDD ~RESET RESET 0 AND Vhigh = 1.8 Vlow = 0 Trise = {tff_rise} Tfall = {tff_rise} Td = td_DFF\r\n'
'.model NMOS NMOS\r\n'...
'.model PMOS PMOS\r\n'...
'.lib C:\\Users\\Gautam Jha\\Documents\\LTspiceXVII\\lib\\cmp\\standard.mos\r\n'...
'.param VDD=1.8V VSS=0V td_DFF=1p T=30n Tp=15n tff_rise=50p\r\n'...
'.ic V(Vp) = 0.9 V(Vpc) = 0.9\r\n'...
'.include tsmc018.txt\r\n'...
'.tran {50*T}\r\n'...
'.LIB opamp.sub\r\n'...
'.param Wn1=4u Ln1=0.25u Wn2=4u Ln2=.25u Wp1=20u Lp1=0.25u Wp2=20u Lp2=0.25u\r\n'...
'.param fn1=8 fn2=8 fp1=8 fp2=8\r\n'...
'.ic V(Va) = 0.9 V(Vb) = 0.9 V(V_SnH) = 0.9\r\n'...
'.backanno\r\n'...
'.end\r\n'];
%}

% fid = fopen('ckt_netlist.txt', 'r');
% %txt_data = textscan(fid, '%s', 'delimiter', '\n');
% txt_data = fscanf(fid, '%s \n');
% fclose(fid);

data = fileread('ckt_netlist2.txt');
txt_data2 = strsplit(data,'\n');
str = string(txt_data2);
string_to_write = ['.param X = 0.4',newline];
str(end - 1:end + 1) = str(end - 2:end);
str(end-3) = string_to_write;

fid = fopen(netlist,'w+');
%fprintf(fid, code);
%fprintf(fid, '%s\r\n ', txt_data);
fprintf(fid, '%s\n', str);
fid = fclose(fid);

fileName = 'Current_Mode_Integrator_3.raw';
dos('LTSpice_call_ckt1.bat');
pause(8) ;
raw_data = LTspice2Matlab(fileName);
dos('LTSpice_end.bat');
% tic

%% Extract sample %%
VDC = 0.9; T = 40e-9; Tp = T/2 ;
Vamp = 0.3; CF = 4; Tp_clk3 = 0.475*T;

id_vin = find(contains(raw_data.variable_name_list,'V(vin)'));
id_vp = find(contains(raw_data.variable_name_list,'V(vp)'));
id_vclk3 = find(contains(raw_data.variable_name_list,'V(vclk3)'));
id_vclk2 = find(contains(raw_data.variable_name_list,'V(vclk2)'));
id_v_snh = find(contains(raw_data.variable_name_list,'V(v_snh)'));

vclk3 = raw_data.variable_mat(id_vclk3, :);
vclk2 = raw_data.variable_mat(id_vclk2, :);
%clk3_edges_id_list = find(vclk3 >= 0.5 & vclk3 <= 0.6);
vp = raw_data.variable_mat(id_vp, :);
%plot(raw_data.time_vect,raw_data.variable_mat(id, :));
vin = raw_data.variable_mat(id_vin, :);
v_snh = raw_data.variable_mat(id_v_snh, :);

len = floor(raw_data.time_vect(end)/(T*CF));
obs_samples_vp = zeros(1, len);

i = (CF - 1)*T + 2*Tp_clk3; j = 1;
while(i < raw_data.time_vect(end))
   sample_id = find(abs(raw_data.time_vect - i) == min(abs(raw_data.time_vect - i))) ;
    obs_samples_vp(j) = vp(sample_id);
   i = i + CF*T; j = j + 1;
end
% obs_samples_vp = obs_samples_vp - VDC;


%% Theta construction %%
N = round(raw_data.time_vect(end)/T);
M = N/CF;


G = 1/(10e3); C = 5e-12;
freq_PNS = 1/T ;
clk2_rise_time = 50e-12;
clk3_rise_time = 50e-12;
clk2_on_time = 0.375*T;


sample_times_clk2 = clk2_rise_time + clk2_on_time : T : raw_data.time_vect(end);
PRBS_mat = construct_PRBS_mat(N, CF);
theta = (Vamp*G*Tp/C)*construct_theta(PRBS_mat, freq_PNS,sample_times_clk2, VDC);
struc_data = load('simdata_N_80.mat','theta_data');
theta_data = struc_data.theta_data - repmat(PRBS_mat * (G*VDC*(ones(N, 1))), [1, N]);
theta_data_id = zeros(size(theta_data));
norm_vals = zeros(1, N);
for i = 1:N
    norm_vals(i) = norm(theta_data(:, i));
    theta_data_id(:, i) = theta_data_id(:, i)/norm_vals(i);
end    
corr_mat_theta = correl(theta);
% theta_data = construct_theta_data(N, CF, fileName, txt_data2, netlist, Tp_clk3, T, VDC, id_vp);
corr_mat_theta_data = correl(theta_data);


obs_samples_vin = zeros(size(sample_times_clk2)); j = 1;
for i = 1:length(sample_times_clk2)
    sample_id = find(abs(raw_data.time_vect - sample_times_clk2(i)) == min(abs(raw_data.time_vect - sample_times_clk2(i)))) ;
    obs_samples_vin(j) = vin(sample_id);
    j = j + 1;
end
obs_samples_vin = obs_samples_vin - VDC;
test_vin = PRBS_mat * obs_samples_vin';


sample_times_clk3 = clk3_rise_time + 2*Tp_clk3 : T : raw_data.time_vect(end);
obs_samples_v_snh = zeros(size(sample_times_clk3)); j = 1;
for i = 1 : length(sample_times_clk3)
    sample_id = find(abs(raw_data.time_vect - sample_times_clk3(i)) == min(abs(raw_data.time_vect - sample_times_clk3(i)))) ;
    obs_samples_v_snh(j) = v_snh(sample_id);
    j = j + 1;
end
obs_samples_v_snh = obs_samples_v_snh - VDC;
test_v_snh = PRBS_mat * obs_samples_v_snh';
error1 = obs_samples_v_snh - obs_samples_vin;


obs_samples_vp = obs_samples_vp - VDC - (PRBS_mat * (G*VDC*(ones(N, 1))))';


%% Test Signal %%
f = 5e6;
test_signal = (0.3*Vamp*G*Tp/C)*sin(2*pi*f*sample_times_clk2);
test_out = PRBS_mat * test_signal';

%% Signal Reconstruction %%

f = linspace(0,freq_PNS/2 - freq_PNS/N, N/2);
% % [sparse_indx,coeff] = dsp_3(theta,obs_samples_vp',1);
% % [sparse_indx,coeff] = dsp_3(theta,test_out,1);
% [sparse_indx,coeff] = dsp_3(theta,(G*Tp/C)*test_vin,1);
% [sparse_indx,coeff] = dsp_3(theta,(G*Tp/C)*test_v_snh,1);
[sparse_indx,coeff] = dsp_4(theta_data_id,theta_data ,obs_samples_vp',1);
sparse_indx
cc = coeff;
sparse_indx = (sparse_indx-1)*freq_PNS/N ;
for i = 2:(N/2)
   coeff(i) = coeff(i) - coeff(N+2-i) ;
end
Reconstructed_spectrum = abs(coeff(1:N/2));

plot(f,Reconstructed_spectrum);
xlabel('f');
ylabel('|X(f)|');
title('Reconstructed spectrum');

% [XX,sparse_indx2] = OMP(1,obs_samples_vp,theta_data);
% sparse_indx2 = sparse_indx2 * freq_PNS/N ;
% XX
% XX = XX(1:N/2);
% f = f([2:N/2,1]);
% plot(f,XX);
% title(['Reconstructed signal FFT (fs = ',num2str(freq_PNS),'  M = N/',num2str(N/M),') ']);
% xlabel('f'); ylabel('|X''(f)|');
% 
% subplot(2, 1, 1);
% plot(raw_data.time_vect, vin); hold on; plot(raw_data.time_vect, vclk2);
% N = 128;
% figure;
% f = linspace(-freq_PNS/2, freq_PNS/2 - freq_PNS/N, N);
% subplot(2,1,1);
% plot(f, abs(fftshift(fft(obs_samples_vin,N)/N)));
% subplot(2,1,2);
% plot(f, abs(fftshift(fft(obs_samples_v_snh ,N)/N)));
% toc

