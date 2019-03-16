clear;clc
load('PN_Code_Baker.mat'); % Load true code (no oversampling version)
filename = [1,1];          % Name of xxdata2.csv, e.g., [1,1] for 11data2.csv
filecode = bi2de(filename,'left-msb');
switch filecode
    case 1
        M = csvread('C:\Users\Han\Google Drive\Blocker Tagging\Jan 22nd\01data2.csv');   % read data
    case 2
        M = csvread('C:\Users\Han\Google Drive\Blocker Tagging\Jan 22nd\10data2.csv');   % read data
    case 3
        M = csvread('C:\Users\Han\Google Drive\Blocker Tagging\Jan 22nd\11data2.csv');   % read data
end
%%

data_length = size(M,1);            % Data length
Ts = M(2,1) - M(1,1);               % Sample duration
time_win = data_length * Ts;        % Time window of capture
t = M(:,1);                         % Time scale
if filecode == 3
    y_waveform_SW = M(:,2);             % Sam added Square-Wave for timing
    y_waveform_noalign = M(:,3);                % True waveform after LNA
    
    % process square wave for timing recovery
    [start_idx,~] = min(find(y_waveform_SW>0.245));
    start_idx_mid = start_idx + round(0.25e-6/Ts); % pick center of PN code by shifting 0.25us
    y_waveform = [y_waveform_noalign(start_idx_mid:end);y_waveform_noalign(1:start_idx_mid-1)];
    
    % Plot of square wave and staring time
%     figure
%     plot(t/1e-6,y_waveform_SW);hold on
%     grid on
%     xlabel('time [us]')
%     plot(t(start_idx)/1e-6,y_waveform_SW(start_idx),'o','linewidth',2)
    
else
    y_waveform = M(:,2);                % True waveform after LNA
end

%% PSD watch for diagnosis
% fft_size = 8192;                    % FFT size 
% YDATA = abs(fft(y_waveform(1:fft_size))).^2;
% XDATA = linspace(-1,1,fft_size)*(1/(2*Ts));
% figure
% semilogy(XDATA/1e6, [YDATA(fft_size/2+1:end);YDATA(1:fft_size/2)]);
% hold on
% grid on
% xlabel('Freq [MHz]')
% ylabel('PSD [dB]')
%% Processing part to extract 4MHz signals
if filecode == 3
    watch_win = 100e-6; % Time window in 11 is 100us
else
    watch_win = 20e-6;  % Time window in other files is 20us
end
watch_size = floor(watch_win/Ts);

% ---- Digital Filter to Reject Self-Mixing (After Downsampling) -------
fcutoff = 0.01;                   % Fcutoff (0.9 is heristic value)
AA_taps_DS = 2000;                % Num of Tap (just for RF sim purpose)
bhi_AA_DS = fir1(AA_taps_DS,fcutoff,'low'); % Call function for filter design

carrier_OS = cos(0.0*pi + 2*pi*(4 + -0.0e-3 )*1e6*t(1:watch_size));
DSP_in_H0_u = y_waveform(1:watch_size);

sig_analog_DS = filter(bhi_AA_DS,1,(DSP_in_H0_u.*carrier_OS));
sig_analog_DS_shift = sig_analog_DS(AA_taps_DS/2+1 : end);

% Decimator filter to better process signal at 4MHz
DS_ratio = 100;
sig_DS = sig_analog_DS_shift(1:DS_ratio:end);

fcutoff = 0.005;                   % Fcutoff (0.9 is heristic value)
AA_taps_DS2 = 4000;                % Num of Tap (just for RF sim purpose)
bhi_AA_DS2 = fir1(AA_taps_DS2,fcutoff,'low'); % Call function for filter design

sig_analog_DS2 = filter(bhi_AA_DS2,1,sig_DS);
sig_analog_DS2_shift = sig_analog_DS2(AA_taps_DS2/2+1:end);
%%
% figure
% plot(t(1:length(sig_analogBB_shift))/1e-6, sig_analogBB_shift)
% xlabel('Time [us]')
% ylabel('Mag')
% grid on
%% plot baseband signal and PN codes 
Newt = t(1:DS_ratio:end);
if filecode==3
    t_sel = 1:200:length(sig_analog_DS2_shift); % when Ts = 2.5e-11
else
    t_sel = 1:400:length(sig_analog_DS2_shift); % when Ts = 1.25e-11
end

figure
subplot(211) % PN * BLK terms down-converted at DC
plot(Newt(1:length(sig_analog_DS2_shift))/1e-6, sig_analog_DS2_shift);hold on
plot(Newt(t_sel)/1e-6, sig_analog_DS2_shift(t_sel),'ro');
xlabel('Time [us]')
ylabel('Mag')
grid on
title(['PN \times BLK after Process.'])
if filecode ==3
    xlim([-50,-25])
else
    xlim([-10,8])
end

subplot(212) % True code
if filecode ==2
    stem(-10:0.5:7,Code(212:212+35-1,1));hold on
    xlim([-10,8])
    grid on
    xlabel('Time [us]')
    title('Detected True PN code')
    ylabel('Mag')
elseif filecode ==1
    stem(-10:0.5:7,Code(7:7+35-1,1));hold on
    xlim([-10,8])
    grid on
    xlabel('Time [us]')
    title('Detected True PN code')
    ylabel('Mag')
elseif filecode==3
    code1 = Code(1:1+190-1,1);
    code2 = Code(201:201+190-1,1);
    alpha1 = (code1'*sig_analog_DS2_shift(t_sel))/norm(code1)^2;
    alpha2 = (code2'*sig_analog_DS2_shift(t_sel))/norm(code2)^2;
    stem(-50:0.5:44.5, (code1*alpha1 + code2*alpha2) )
    xlim([-50,25])
    grid on
    xlabel('Time [us]')
    title('Detected True PN \times BLK')
    ylabel('Mag')
end


%%
PN_code = zeros(length(sig_analog_DS2_shift(t_sel)),1);
PN_code(find(sig_analog_DS2_shift(t_sel)>0)) = 1;
PN_code(find(sig_analog_DS2_shift(t_sel)<0)) = -1;

figure
stem(PN_code,'o');hold on
grid on
%% Exhaustive search for code words (due to lack of timing)
if filecode == 3
    PN_start_idx = 1;
    figure;
    subplot(211)
    corr_out1 = abs(cconv(code1,flipud(sig_analog_DS2_shift(t_sel))));hold on
    plot(corr_out1);grid on
    plot(190, corr_out1(190), 'o');grid on
    xlabel('Time Shift')
    ylabel('Mag')
    legend('Corr w/ all possible segment','Corr w/ true segment of PN1')
    title('output of correlation w/ PN code 1')
    subplot(212)
    corr_out2 = abs(cconv(code2,flipud(sig_analog_DS2_shift(t_sel))));hold on
    plot(corr_out2);grid on
    plot(190, corr_out2(190), 'o');grid on
    xlabel('Time Shift')
    ylabel('Mag')
    legend('Corr w/ all possible segment','Corr w/ true segment of PN2')
    title('output of correlation w/ PN code 2')
elseif filecode == 1
%     figure;plot(abs(cconv(Code(:,1),flipud(PN_code))));grid on
    PN_start_idx = 7;
    figure;plot(abs(cconv(Code(PN_start_idx:PN_start_idx+35-1,1),flipud(PN_code))));grid on
elseif filecode == 2
    PN_start_idx = 212;
    figure;plot(abs(cconv(Code(PN_start_idx:PN_start_idx+35-1,1),flipud(PN_code))));grid on
end
%% PSD watch
% YDATA = abs(fft(y_waveform(1:fft_size))).^2;
% XDATA = linspace(-1,1,fft_size)*(1/(2*Ts));
% figure
% semilogy(XDATA/1e6, [YDATA(fft_size/2+1:end);YDATA(1:fft_size/2)]);
% hold on
% grid on
% xlabel('Freq [MHz]')
% ylabel('PSD [dB]')