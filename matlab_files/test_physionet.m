clear;

% if ~isfile('PhysionetData.mat')
%     ReadPhysionetData   
% end
% load PhysionetData
% ecg = Signals{10, 1};
% plot(ecg);

[t,signal,Fs,labels]=rdmat('rec_7m');
ecg_data = signal(:, 2)';
% t = tm;
plot(t, ecg_data);
