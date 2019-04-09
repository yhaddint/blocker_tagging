clear;clc
M1 = addpath('C:\Users\Han\Google Drive\Blocker Tagging\March 27th');   % read data
load('data1.mat');
%%
fc = 900.1;
fcutoff = 0.001;                   % Fcutoff (0.9 is heristic value)
AA_taps_DS = 4000;                % Num of Tap (just for RF sim purpose)
bhi_AA_DS = fir1(AA_taps_DS,fcutoff,'low'); % Call function for filter design

Ts = XDelta;
watch_win = 100e-6;
watch_size = floor(watch_win/Ts);
t = (0:watch_size-1).' * Ts;
Ts = XDelta;
carrier_OS = cos(2*pi*(fc + -0.35e-3 )*1e6*t(1:watch_size));

sig_analog_DS = filter(bhi_AA_DS,1,(Y(1:watch_size).*carrier_OS));
sig_analog_DS_shift = sig_analog_DS(AA_taps_DS/2+1 : end);

%%
start_idx = 1;
t_sel = round(start_idx:(6.25e-9/Ts):length(sig_analog_DS_shift));
figure
plot(t(1:length(sig_analog_DS_shift))/1e-6, real(sig_analog_DS_shift));hold on
plot(t(t_sel)/1e-6, real(sig_analog_DS_shift(t_sel)),'o','linewidth',2)
xlabel('Time [us]')
title('PN sequence from data1.mat')
% xlim([0,20])
grid on
%%
seq = sig_analog_DS_shift(t_sel);
figure
plot(-99.56 + t(t_sel(1:1.5e4))/1e-6,cos(1.33 + 2*pi*100e3 * t(t_sel(1:1.5e4))).*real(seq(1:1.5e4)))
xlabel('time [us]')
grid on
%%
filename = [1,0];          % Name of xxdata2.csv, e.g., [1,1] for 11data2.csv
filecode = bi2de(filename,'left-msb');
switch filecode
    case 1
        M = csvread('C:\Users\Han\Google Drive\Blocker Tagging\March 27th\01_10M.csv');   % read data
    case 2
        M = csvread('C:\Users\Han\Google Drive\Blocker Tagging\March 27th\10_10M.csv');   % read data
    case 3
        M = csvread('C:\Users\Han\Google Drive\Blocker Tagging\March 27th\11_10M.csv');   % read data
end
%

data_length = size(M,1);            % Data length
Ts = M(2,1) - M(1,1);               % Sample duration
time_win = data_length * Ts;        % Time window of capture
t_new = M(:,1);                         % Time scale
y_waveform = M(:,2)-mean(M(:,2));
%
figure
plot(M(1:1.5e4,1)/1e-6,y_waveform(1:1.5e4));grid on
%%
freq_cand_num = 40;
freq_cand = linspace(-8,8,freq_cand_num);
phi_cand_num = 100;
phi_cand = linspace(0,pi,phi_cand_num);
% for ff=1:freq_cand_num
%     dfreq = freq_cand(ff)*1e3;
    for ii=1:phi_cand_num
        phi = phi_cand(ii);
        wave_cand = cos(phi + 2*pi * (100e3 + -2.31e3) * t(t_sel(1:1.5e4))).*real(seq(1:1.5e4));
        score(ii) = abs(sum(wave_cand.*y_waveform(1:1.5e4)));
    end
%     score_ff(ff) = max(score);
% end
figure
plot(phi_cand,score);grid on
%%
wave_cand = cos(0 + 2*pi*(100e3-0e3) * t(t_sel(1:1.5e4))).*real(seq(1:1.5e4));
alpha = pinv(wave_cand) * y_waveform(1:1.5e4);
figure
plot(-99.56 + t(t_sel(1:1.5e4))/1e-6, 1.2*alpha * wave_cand,'linewidth',2);
hold on
plot(-99.56 + t(t_sel(1:1.5e4))/1e-6, y_waveform(1:1.5e4),'linewidth',2)
set(gca,'FontSize',14)
grid on
legend('True PN \times cos(\omega t+\phi)','Actual Capture (Subtract DC)')
xlabel('Time [us]')
ylabel('Mag')
% xlim([-75,-65])

%%
car_hat = cos(0 + 2*pi*(100e3-0e3) * t(t_sel(1:1.5e4)));
waveform_true = y_waveform(1:1.5e4).*(car_hat./(abs(car_hat).^2+1e-3));
figure
plot(-99.56 + t(t_sel(1:1.5e4))/1e-6, waveform_true/alpha);hold on
plot(-99.56 + t(t_sel(1:1.5e4))/1e-6, real(seq(1:1.5e4)),'linewidth',2)
% ylim([-1,1])
xlabel('Time [us]')
ylabel('Mag')
legend('Recovered PN from capture','PN Ground Truth')
grid on
% xlim([-75,-65])
%% Correlation in the IF band
scores = abs(sum((wave_cand) .* (y_waveform(1:1.5e4))))
