function [theta_data] = construct_theta_data(N, CF, fileName, txt_data2, netlist, Tp_clk3, T, VDC, id_vp)

M = N/CF;
theta_data = zeros(M, N);
val = 0 : 1/N : (1 - 1/N);
obs_samples_vp = zeros(1, M);

for i = 1:length(val)
    str = string(txt_data2);
    string_to_write = ['.param X =',num2str(val(i)),newline];
    str(end - 1:end + 1) = str(end - 2:end);
    str(end-3) = string_to_write;
    
    fid = fopen(netlist,'w+');
    fprintf(fid, '%s\n', str);
    fid = fclose(fid);
    dos('LTSpice_call_ckt1.bat');
    pause(12) ;
    raw_data = LTspice2Matlab(fileName);
    dos('LTSpice_end.bat');
    
    vp = raw_data.variable_mat(id_vp, :);
    k = (CF - 1)*T + 2*Tp_clk3; j = 1;
    while(k < raw_data.time_vect(end))
      sample_id = find(abs(raw_data.time_vect - k) == min(abs(raw_data.time_vect - k))) ;
      obs_samples_vp(j) = vp(sample_id(1));
      k = k + CF*T; j = j + 1;
    end
    obs_samples_vp = obs_samples_vp - VDC;
    
    theta_data(:, i) = obs_samples_vp';
end

end

