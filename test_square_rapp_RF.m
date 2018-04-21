% The Rapp approximation model should work in baseband. But it is not clear
% whether it can be used to simulate RF-self mixing of square-law devices

% ------ script control parameter --------
clear;clc;
rng(2)
plot_PSD = 1;
plot_rapp = 1;

% ------ system parameters ---------
sig_length = 1e5;
CAL_num = 2;
CAL_pow = 0.02;
PN1 = randi(2,5e1,1)*2-3;
BLK_BB1 = 1*ones(sig_length,1);
CAL_BB1 = 0.1*kron(PN1,ones(2e3,1));
sig_BB1 = BLK_BB1 + CAL_BB1;
sig1 = sqrt(CAL_pow*2) * cos(pi*2*(0.25/10)*(0:sig_length-1).').*sig_BB1;

PN2 = randi(2,5e1,1)*2-3;
BLK_BB2 = 1*ones(sig_length,1);
CAL_BB2 = 0.1*kron(PN2,ones(2e3,1));
sig_BB2 = BLK_BB2 + CAL_BB2;
sig2 = sqrt(CAL_pow*2) * cos(pi*2*(0.2/10)*(0:sig_length-1).').*sig_BB2;
%% plot input-output relationship of square-law device
if plot_rapp
    figure(100)
    P_range = [0.75,1,2,10];
    xin = 10.^(linspace(-20,5,1e3)/10);
    for pp=1:length(P_range)
        xout = get_rapp_square(xin,1,P_range(pp));
        loglog(xin,xout,'linewidth',2);hold on
    end
    loglog(xin,xin.^2,'k--','linewidth',2);hold on
    grid on
    xlabel('Input Magnitude')
    ylabel('Output Magnitude')
    legend('P = 0.75','P = 1', 'P = 2', 'P = 10','Perfect Square')
    ylim([0,1.1])
end
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
temp_shift = temp(2001:end)-CAL_pow*CAL_num;
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
plot(DSP_in);grid on
subplot(212)
plot(PN1)
%% Correlation & LPF
bhi_ED = fir1(40,0.05,'low');
corr_out = DSP_in.*PN1(1:end-1);
ED_out_temp = filter(bhi_ED,1,corr_out);
ED_out = sum(abs(ED_out_temp(21:end)).^2)
%% time domain signal plot
% figure
% plot(real(sig));hold on
% plot(real(sig_out))