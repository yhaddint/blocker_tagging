clear;clc
M = csvread('data/Time_50us.csv');   % read data
%%

data_length = size(M,1);            % Data length
Ts = M(2,1) - M(1,1);               % Sample duration
time_win = data_length * Ts;        % Time window of capture
t = M(:,1);                         % Time scale
y_waveform = M(:,2);                % waveform after LNA
fft_size = 8192;                    % FFT size 
%% PSD watch
% YDATA = abs(fft(y_waveform(1:fft_size))).^2;
% XDATA = linspace(-1,1,fft_size)*(1/(2*Ts));
% figure
% semilogy(XDATA/1e6, [YDATA(fft_size/2+1:end);YDATA(1:fft_size/2)]);
% hold on
% grid on
% xlabel('Freq [MHz]')
% ylabel('PSD [dB]')
%%
watch_win = 80e-6;
watch_size = floor(watch_win/Ts);
% ---- Digital Filter to Reject Self-Mixing (After Downsampling) -------
fcutoff = 0.01;                   % Fcutoff (0.9 is heristic value)
AA_taps_DS = 2000;                % Num of Tap (just for RF sim purpose)
bhi_AA_DS = fir1(AA_taps_DS,fcutoff,'low'); % Call function for filter design

carrier_OS = cos(2*pi*(210 + -0.45e-3 )*1e6*t(1:watch_size));
DSP_in_H0_u = y_waveform(1:watch_size);

sig_analogBB = filter(bhi_AA_DS,1,(DSP_in_H0_u.*carrier_OS));
sig_analogBB_shift = sig_analogBB(AA_taps_DS/2+1 : end);

figure
plot(t(1:length(sig_analogBB_shift))/1e-6, sig_analogBB_shift)
xlabel('Time [us]')
ylabel('Mag')
grid on

%% PSD watch
YDATA = abs(fft(y_waveform(1:fft_size))).^2;
XDATA = linspace(-1,1,fft_size)*(1/(2*Ts));
figure
semilogy(XDATA/1e6, [YDATA(fft_size/2+1:end);YDATA(1:fft_size/2)]);
hold on
grid on
xlabel('Freq [MHz]')
ylabel('PSD [dB]')