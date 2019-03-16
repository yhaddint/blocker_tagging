clear;clc
M = addpath('C:\Users\Han\Google Drive\Blocker Tagging\Jan 22nd');   % read data
load('data1.mat');
%%
fc = 510;
fcutoff = 0.01;                   % Fcutoff (0.9 is heristic value)
AA_taps_DS = 2000;                % Num of Tap (just for RF sim purpose)
bhi_AA_DS = fir1(AA_taps_DS,fcutoff,'low'); % Call function for filter design

Ts = XDelta;
watch_win = 100e-6;
watch_size = floor(watch_win/Ts);
t = (0:watch_size-1).' * Ts;
Ts = XDelta;
carrier_OS = cos(2*pi*(fc + -0.45e-3 )*1e6*t(1:watch_size));

sig_analog_DS = filter(bhi_AA_DS,1,(Y(1:watch_size).*carrier_OS));
sig_analog_DS_shift = sig_analog_DS(AA_taps_DS/2+1 : end);

%%
t_sel = 1e4:round(0.5e-6/Ts):length(sig_analog_DS_shift);
figure
plot(t(1:length(sig_analog_DS_shift)), real(sig_analog_DS_shift));hold on
plot(t(t_sel), real(sig_analog_DS_shift(t_sel)),'o','linewidth',2)
grid on
%%
PN_code = zeros(length(real(sig_analog_DS_shift(t_sel))),1);
PN_code(find(real(sig_analog_DS_shift(t_sel))>0)) = 1;
PN_code(find(real(sig_analog_DS_shift(t_sel))<0)) = -1;
% figure;stem(PN_code)

load('PN_Code_Baker.mat');
figure;plot(abs(cconv(Code(:,1),flipud(PN_code))));grid on

%%
if fc==210
    start_idx = 1;
    sum(PN_code.*Code(start_idx:start_idx+length(PN_code)-1,1))
    figure;plot(abs(conv(Code(start_idx:start_idx+200-1,1),flipud(PN_code))));grid on
elseif fc==510
    start_idx = 201;
    sum(PN_code.*Code(start_idx:start_idx+length(PN_code)-1,1))
    figure;plot(abs(conv(Code(start_idx:start_idx+200-1,1),flipud(PN_code))));grid on
end
%%
% fft_size = 8192;                    % FFT size 
% YDATA = abs(fft(Y(1:fft_size))).^2;
% XDATA = linspace(-1,1,fft_size)*(1/(2*Ts));
% figure
% semilogy(XDATA/1e9, [YDATA(fft_size/2+1:end);YDATA(1:fft_size/2)]);
% hold on
% grid on
% xlabel('Freq [GHz]')
% ylabel('PSD [dB]')