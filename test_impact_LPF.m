clear;clc;
seq = randi(2,100,1)*2-3;
seq_intp = kron(seq,ones(10,1));
%%
% figure
% plot(seq_intp)
XDATA = linspace(-1,1,1e3).';
YDATA = 20*log10(abs(fft(seq_intp)));
figure
plot(XDATA, [YDATA(501:1000);YDATA(1:500)])
grid on
xlabel('Freq [\pi]')

%% LPF
fcutoff = 0.4;                   % Fcutoff (200MHz is heristic value)
AA_taps_LP = 10;                % Num of Tap (just for RF sim purpose)
bhi_AA_LP = fir1(AA_taps_LP, fcutoff, 'low'); % Call function for filter design

%%
seq_filter = filter(bhi_AA_LP,1,[seq_intp;zeros(5,1)]);
carrier1 = cos(0.1/10*2*pi*(1:1e3)).';
carrier1_hat = sin(0.00*pi + 0.1/10*2*pi*(1:1e3)).';

% XDATA = linspace(-1,1,1e3).';
YDATA = 20*log10(abs(fft(carrier1 .* carrier1_hat .* seq_filter(6:end))));
figure
plot(XDATA, [YDATA(501:1000);YDATA(1:500)])
grid on
xlabel('Freq [\pi]')
%%
% for ii=1:10
%     phi = 2*pi/10*(ii-1+0.02);
%     H = diag(cos(phi + 0.1/10*2*pi*(1:1e3).'));
%     results(:,ii) = pinv(H)* (carrier1 .* seq_filter(6:end));
% end
% %%
% figure
% plot(results(:,1))
% 
%%
figure
plot(carrier1.*seq_filter(6:end));
grid on

