% The Rapp approximation model should work in baseband. But it is not clear
% whether it can be used to simulate RF-self mixing of square-law devices

% ------ script control parameter --------
clear;clc;
rng(2)
plot_PSD = 1;

% ------ system parameters ---------
sig_length = 1e5;
PN1 = randi(2,5e1,1)*2-3;
sig_BB1 = kron(PN1+1,ones(2e3,1));
sig1 = 0.1 * cos(pi*2*(0.25/10)*(0:sig_length-1).').*sig_BB1;

PN2 = randi(2,5e1,1)*2-3;
sig_BB2 = kron(PN2+1,ones(2e3,1));
sig2 = 0.1 * cos(pi*2*(0.2/10)*(0:sig_length-1).').*sig_BB2;

%% Rapp model for square-law devices
sig = sig1+sig2;
sig_out = get_rapp_square(sig,1,1);


%% frequency domain signal plot
if plot_PSD
    figure
    [ x_freq, PSD_actual ] = get_PSD( sig, 4*8192, 1e4 );
    subplot(311)
    plot(x_freq, PSD_actual);hold on
    xlim([0,1e3])
    grid on
    xlabel('Freq (MHz)')
    ylabel('PSD')
    [ x_freq, PSD_actual ] = get_PSD( sig_out, 4*8192, 1e4 );
    subplot(312)
    plot(x_freq, PSD_actual);hold on
    xlim([0,1e3])
    grid on
    xlabel('Freq (MHz)')
    ylabel('PSD')

    subplot(313)
    plot(x_freq, PSD_actual);hold on
    xlim([-1e1,1e1])
    grid on
    xlabel('Freq (MHz)')
    ylabel('PSD')
end
%% Baseband signal (after anti-aliasing filter but before ADC)

bhi = fir1(4000,0.003,'low');
temp = filter(bhi,1,sig_out);
temp_shift = temp(2001:end);
xdata = linspace(0,1e5/10e9*1e6,length(temp_shift));
figure(99)
plot(xdata,temp_shift);grid on
xlabel('Time (\mu s)')
ylabel('Magnitude')

%% Sampling with aux ADC (~10MHz)
DSP_in = downsample(temp_shift(1e3:end),2e3);
DSP_time_index = downsample(xdata(1e3:end),2e3);
figure(98)
subplot(211)
plot(DSP_time_index,DSP_in-mean(DSP_in));grid on
subplot(212)
plot(PN2)

%% time domain signal plot
% figure
% plot(real(sig));hold on
% plot(real(sig_out))