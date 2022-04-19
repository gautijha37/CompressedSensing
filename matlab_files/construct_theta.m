function [theta] = construct_theta(PRBS_mat, freq_PNS, sample_times, VDC)
    N = length(sample_times);
    f = 0 : freq_PNS/N : (N - 1)/N * freq_PNS ;
    x = sin(2*pi* sample_times' * f);
    theta = PRBS_mat*x;
end

