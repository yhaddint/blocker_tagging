

% ------ script control parameter --------
clear;clc;warning off
rng(2)                              % Random seed
plot_PSD = 1;                       % Whether plot PSD figure
plot_rapp = 0;                      % Whether plot IP/OP relationship in LNA 
MCtimes = 1e0;                    % Number of Monte Carlo trials

% ------ system parameters ---------
T_win = 500e-6;                     % Time window to capture signal in [s], i.e., 200us 
fs = 5;                             % RF Nyquist rate for simulator [GHz]
PN_BW = 0.002;                      % BW of PN signal in [GHz], i.e.,2MHz
PN_OS = fs/PN_BW;                   % Oversampling ratio for PN signal
BLK_BW = 0.0001;                    % BW of blocker in [GHz], i.e., 100KHz
BLK_OS = fs/BLK_BW;                 % Oversampling ratio for blocker signal
sig_length = T_win*fs*1e9;          % Number of sample (conceptual sample in RF domain)
BLK_num = 2;
% CAL_pow = 0.04;                   % Power of calibration (PN) signal (to model saturation)
BLK_pow = [0.1 0.1];                % Blocker power (relative to PN power)
freq = [0.2,0.5];                   % Center frequency of blockers [GHz]
mag_H1 = [1, 1];                    % Enable/disable blockers in hypothesis testing
mag_H0 = [0, 1];                    % Enable/disable blockers in hypothesis testing
resistance = 50;                    % default resistance [om]
PN_CFO = 10e6/1e9;                  % freq. offset b/w each PN and blocker in [GHz]

% setting when PN and BLK has freq offset
if PN_CFO~=0
    ADC_OS = 5;                     % ADC need to be faster the such offset; Default 5x faster than freq_OS
    shifted_PN_OS = fs/PN_CFO;      % Tagging Aux. ADC sampling rate in terms of fs
end


% ---- Sample Sequences of Blockers -------
% Sim. uses 2 stages - 1) 50x oversampling with pulse shaping function
%                      2) 500x interpolation without pulse shaping
MQAM = 16;                              % Modulation level, 16QAM
upsam = 50;                             % Oversampling ratio in stage 1
upsam_addition = BLK_OS/upsam;          % Requires interpolation in stage 2
BLK_symb_temp = real( get_SC_waveform( MQAM, 1e5, upsam, BLK_num ) );
BLK_symb = BLK_symb_temp./norm(BLK_symb_temp,'fro')*sqrt(sig_length/10*BLK_num);

% ---- Analog Anti-aliasing Filter -------
fcutoff = 0.040;                          % Fcutoff freq in [GHz], i.e., 40 MHz. Greater than (PN BW + PN/BLK freq offset)
fcutoff_DSP = fcutoff/fs*2;               % Effective "digital" frequency in RF sim.
AA_taps = 2000;                           % Num of Tap (just for RF sim purpose)
bhi_AA = fir1(AA_taps,fcutoff_DSP,'low'); % Call function for filter design

% ---- Digital Filter to Reject Self-Mixing (After Downsampling) -------
fcutoff = 2/ADC_OS*0.9;                   % Fcutoff (0.9 is heristic value)
AA_taps_DS = 100;                         % Num of Tap (just for RF sim purpose)
bhi_AA_DS = fir1(AA_taps_DS,fcutoff,'low'); % Call function for filter design

% ---- Filter to reject PN ------

LPF_cut_margin = BLK_BW/PN_BW*1.25;       % LPF in DSP (keep blocker & reject PN spreaded blockers, 1.25 is heuristic value)
PN_reject_taps = 40;
bhi_LPF = fir1(PN_reject_taps,LPF_cut_margin,'low');  % Filter design

% pow_range = [2.5,0,-2.5,-5,-7.5,-10,-15,-25,-35]; % Test various IP power
pow_range = [-20];                                % Total input power range in [dBm] (BLK and PN)

% Zero initialization of matrces
ED_out_H0 = zeros(length(pow_range),MCtimes);
ED_out_H1 = zeros(length(pow_range),MCtimes);

% FOR LOOP for input power
for pp=1:length(pow_range)
    
    IP_mag = dBm2V( pow_range(pp),resistance ); % Input magnitude in [Volt]
    
    % FOR LOOP for Monte Carlo realizations
    for MCindex = 1:MCtimes
       
        % FOR LOOP of each band in a multiband signal
        for bb=1:BLK_num
            
            % --- Pick signal from a T_win time window ------ 
            start_index = randi(length(BLK_symb)-sig_length/upsam_addition)-1;
            BLK_capture(:,bb) = BLK_symb(start_index+1:start_index+sig_length/upsam_addition,bb);
            time_orig = linspace(0,sig_length-1,sig_length/upsam_addition);
            time_intep = linspace(0,sig_length-1,sig_length);
            BLK_BB(:,bb) = sqrt(BLK_pow(bb))*interp1(time_orig, BLK_capture(:,bb), time_intep);

            % ---- PN sequence -------
            % Interpolation of \pm1 sequences
            PN_u(:,bb) = randi(2,sig_length/PN_OS,1)*2-3;
            PN_l(:,bb) = randi(2,sig_length/PN_OS,1)*2-3;
            CAL_BB_u(:,bb) = kron(PN_u(:,bb),ones(PN_OS,1));
            CAL_BB_l(:,bb) = kron(PN_l(:,bb),ones(PN_OS,1));

            % ---- In-band signal (baseband version) -------
            sig_BB_H0_u(:,bb) = mag_H0(bb) * BLK_BB(:,bb);% + CAL_BB_u(:,bb);
            sig_BB_H0_l(:,bb) = mag_H0(bb) * BLK_BB(:,bb);% + CAL_BB_l(:,bb);
            sig_BB_H1_u(:,bb) = mag_H1(bb) * BLK_BB(:,bb);% + CAL_BB_u(:,bb);
            sig_BB_H1_l(:,bb) = mag_H1(bb) * BLK_BB(:,bb);% + CAL_BB_l(:,bb);
            
            % Carrier (offset b/w PN and BLK)
            carrier_BLK(:,bb) = cos(pi*2*(freq(bb)/fs)*(0:sig_length-1).');
            carrier_PN(:,bb) = cos(pi*2*((freq(bb)+PN_CFO)/fs)*(0:sig_length-1).');

            % ---- RF signals ------
            sig_H0_u(:,bb) = IP_mag * (carrier_BLK(:,bb) .* sig_BB_H0_u(:,bb)...
                                      +carrier_PN(:,bb) .* CAL_BB_u(:,bb));
                                  
            sig_H0_l(:,bb) = IP_mag * (carrier_BLK(:,bb) .* sig_BB_H0_l(:,bb)...
                                      +carrier_PN(:,bb) .* CAL_BB_l(:,bb));
                                  
            sig_H1_u(:,bb) = IP_mag * (carrier_BLK(:,bb) .* sig_BB_H1_u(:,bb)...;
                                     +carrier_PN(:,bb) .* CAL_BB_u(:,bb));
            
            sig_H1_l(:,bb) = IP_mag * (carrier_BLK(:,bb) .* sig_BB_H1_l(:,bb)...
                                      +carrier_PN(:,bb) .* CAL_BB_l(:,bb));
            
        end
        
        % ---- Input RF Waveform of LNA -------
        sig_in_H0_u = sum(sig_H0_u,2);
        sig_in_H0_l = sum(sig_H0_l,2);
        sig_in_H1_u = sum(sig_H1_u,2);
        sig_in_H1_l = sum(sig_H1_l,2);
        IP_RMS_pow_H0 = V2dBm( sqrt(mean(abs(sig_in_H0_u).^2)) );
        IP_RMS_pow_H1 = V2dBm( sqrt(mean(abs(sig_in_H1_u).^2)) );
        
        % --- print iteration number ----
        clc;fprintf('IP Power Ite %d outof %d\nMC Ite %d out of %d\n',...
                      pp, length(pow_range), MCindex, MCtimes);
        fprintf('IP power (requested) %4.2f dBm\n', pow_range(pp));
        fprintf('IP power (actual in H0) %4.2f dBm\n', IP_RMS_pow_H0);
        fprintf('IP power (actual in H1) %4.2f dBm\n', IP_RMS_pow_H1);

        % ----- Simulation of square-law devices in RF domain -----
        Asat = 1;
        P_rapp = 1;
        DC_u = 0;%(rand*2-1)*1e-6;
        DC_l = 0;%(rand*2-1)*1e-6;
        sig_out_H0_u = get_rapp_square( sig_in_H0_u, Asat, P_rapp ) + DC_u;
        sig_out_H0_l = get_rapp_square( sig_in_H0_l, Asat, P_rapp ) + DC_l;
        sig_out_H1_u = get_rapp_square( sig_in_H1_u, Asat, P_rapp ) + DC_u;
        sig_out_H1_l = get_rapp_square( sig_in_H1_l, Asat, P_rapp ) + DC_l;

        % ---- Anti-Aliasing LPF before aux ADC -----
        AA_shift = AA_taps/2+1; % Delay after AA filter
        temp_H0_u = filter(bhi_AA,1,sig_out_H0_u);
        temp_H0_l = filter(bhi_AA,1,sig_out_H0_l);
        temp_H1_u = filter(bhi_AA,1,sig_out_H1_u);
        temp_H1_l = filter(bhi_AA,1,sig_out_H1_l);
        
        
        % signal processing starts here      
        if PN_CFO==0 % Original setting (PN on top of BLK)
            
            % ----- Compensate signals from PN self-mixing in DC --------
            sig_analogBB_H0_u = (temp_H0_u(AA_shift:end)...
                                    -(IP_mag^2)/2*(BLK_num))*sqrt(1/IP_mag^2);
            sig_analogBB_H0_l = (temp_H0_l(AA_shift:end)...
                                    -(IP_mag^2)/2*(BLK_num))*sqrt(1/IP_mag^2);
            sig_analogBB_H1_u = (temp_H1_u(AA_shift:end)...
                                    -(IP_mag^2)/2*(BLK_num))*sqrt(1/IP_mag^2);
            sig_analogBB_H1_l = (temp_H1_l(AA_shift:end)...
                                    -(IP_mag^2)/2*(BLK_num))*sqrt(1/IP_mag^2);
            
            % ------ Sampling with aux ADC (~same rate as baseband PN) --------
            corr_in_H0_u = downsample( sig_analogBB_H0_u(PN_OS/2:end), PN_OS );
            corr_in_H0_l = downsample( sig_analogBB_H0_l(PN_OS/2:end), PN_OS );
            corr_in_H1_u = downsample( sig_analogBB_H1_u(PN_OS/2:end), PN_OS );
            corr_in_H1_l = downsample( sig_analogBB_H1_l(PN_OS/2:end), PN_OS );
            
        else % Sam's setting (PN near BLK by freq_OS)
            
            % ------ Sampling with aux ADC (~5X rate as PN/BLK offset) --------
            DSP_in_H0_u = downsample( temp_H0_u(AA_shift:end), shifted_PN_OS/ADC_OS );
            DSP_in_H0_l = downsample( temp_H0_l(AA_shift:end), shifted_PN_OS/ADC_OS );
            DSP_in_H1_u = downsample( temp_H1_u(AA_shift:end), shifted_PN_OS/ADC_OS );
            DSP_in_H1_l = downsample( temp_H1_l(AA_shift:end), shifted_PN_OS/ADC_OS );
            
            % ---- Shift (PN times BLK) to DC -------
            carrier_OS = cos(pi*2*(1/ADC_OS)*(0:length(DSP_in_H0_u)-1).');
            sig_analogBB_H0_u_shift = filter(bhi_AA_DS,1,(DSP_in_H0_u.*carrier_OS));
            sig_analogBB_H0_l_shift = filter(bhi_AA_DS,1,(DSP_in_H0_l.*carrier_OS));
            sig_analogBB_H1_u_shift = filter(bhi_AA_DS,1,(DSP_in_H1_u.*carrier_OS));
            sig_analogBB_H1_l_shift = filter(bhi_AA_DS,1,(DSP_in_H1_l.*carrier_OS));
            
            AA_shift_DS = AA_taps_DS/2+1;
            sig_analogBB_H0_u = sig_analogBB_H0_u_shift(AA_shift_DS:end);
            sig_analogBB_H0_l = sig_analogBB_H0_l_shift(AA_shift_DS:end);
            sig_analogBB_H1_u = sig_analogBB_H1_u_shift(AA_shift_DS:end);
            sig_analogBB_H1_l = sig_analogBB_H1_l_shift(AA_shift_DS:end);
            
            % --------- Down Sample to PN rate -------------
            corr_in_H0_u = downsample( sig_analogBB_H0_u, ADC_OS*PN_CFO/PN_BW );
            corr_in_H0_l = downsample( sig_analogBB_H0_l, ADC_OS*PN_CFO/PN_BW );
            corr_in_H1_u = downsample( sig_analogBB_H1_u, ADC_OS*PN_CFO/PN_BW );
            corr_in_H1_l = downsample( sig_analogBB_H1_l, ADC_OS*PN_CFO/PN_BW );
            
        end
        
        % ---- PN Despreading (Assuming Perfect Timing b/w LNA IP/OP) ----
        corr_out_H0_u = corr_in_H0_u.*PN_u(1:length(corr_in_H0_u),1);
        corr_out_H0_l = corr_in_H0_l.*PN_l(1:length(corr_in_H0_l),1);
        corr_out_H1_u = corr_in_H1_u.*PN_u(1:length(corr_in_H1_u),1);
        corr_out_H1_l = corr_in_H1_l.*PN_l(1:length(corr_in_H1_l),1);

        % ---- LPF to Reject Spreaded Signals ----
        PN_reject_shift = PN_reject_taps/2+1;
        
        LPF_out_temp_H0_u = filter(bhi_LPF,1,corr_out_H0_u);
        LPF_out_temp_H0_l = filter(bhi_LPF,1,corr_out_H0_l);
        LPF_out_temp_H1_u = filter(bhi_LPF,1,corr_out_H1_u);
        LPF_out_temp_H1_l = filter(bhi_LPF,1,corr_out_H1_l);
        
        LPF_out_H0_u = LPF_out_temp_H0_u(PN_reject_shift:end);
        LPF_out_H0_l = LPF_out_temp_H0_l(PN_reject_shift:end);
        LPF_out_H1_u = LPF_out_temp_H1_u(PN_reject_shift:end);
        LPF_out_H1_l = LPF_out_temp_H1_l(PN_reject_shift:end);

        % ----- Cross-Multiplication and ED -------
        ED_out_H0(MCindex,pp) = mean(LPF_out_H0_u.*LPF_out_H0_l);
        ED_out_H1(MCindex,pp) = mean(LPF_out_H1_u.*LPF_out_H1_l);
       
    end % end of Monte Carlo loop
end % end of IP power sweeping

%% Evaluation the detection statistic (for debug)
% pp=1;
% figure(99)
% [b,a] = ksdensity(ED_out_H0(:,pp));
% plot(a,b,'linewidth',2);hold on
% [b,a] = ksdensity(ED_out_H1(:,pp));
% plot(a,b,'linewidth',2);hold on
% grid on
% xlabel('ED output');
% ylabel('PDF');
% legend('H0','H1')
% % detection statistic in H0; Based on 2016/09 slides page 20
% CAL_pow(pp) = 10^((pow_range(pp)-30)/10)*100;
% sigma2_n = LPF_cut_margin * (CAL_pow(pp)*BLK_pow(2)); 
% mu_H0 = 0;
% N_effect = sig_length/BLK_OS;
% sigma2_H0 = sigma2_n^2/N_effect;
% x = linspace(mu_H0 - 4*sqrt(sigma2_H0), mu_H0 + 4*sqrt(sigma2_H0), 1e3);
% y = normpdf(x,mu_H0,sqrt(sigma2_H0));
% figure(99)
% plot(x,y,'--','linewidth',2);hold on
% % detection statistic in H; Based on 2016/09 slides page 20
% sigma2_s = CAL_pow(pp)*BLK_pow(1)*0.85; 
% mu_H1 = sigma2_s;
% N_effect = sig_length/BLK_OS;
% sigma2_H1 = (2*sigma2_s^2 + 2*sigma2_s*sigma2_n + sigma2_n^2)/N_effect;
% % sigma2_H1 = (2*sigma2_s^2)/N_effect;
% 
% x = linspace(mu_H1 - 4*sqrt(sigma2_H1), mu_H1 + 4*sqrt(sigma2_H1), 1e3);
% y = normpdf(x,mu_H1,sqrt(sigma2_H1));
% figure(99)
% plot(x,y,'--','linewidth',2);hold on
%% PSD evalution (for debug)
% plot_PSD = 1;
% if plot_PSD
%     figure
%     [ x_freq, PSD_actual ] = get_PSD( corr_out_H1_u, 998, 40 );
%     subplot(211)
%     plot(x_freq, PSD_actual);hold on
%     xlim([-20,20])
%     grid on
%     xlabel('Freq.')
%     ylabel('PSD')
%     [ x_freq, PSD_actual ] = get_PSD( corr_out_H1_l, 998, 40 );
%     subplot(212)
%     plot(x_freq, PSD_actual);hold on
%     xlim([-20,20])
%     grid on
%     xlabel('Freq.')
%     ylabel('PSD')
% 
% end
%%  Detection performance evaluation 
if MCtimes>1
    if PN_CFO == 0 % uses the theoretical optimal TH
        for pp = 1:length(pow_range)
            CAL_pow = 10^((pow_range(pp)-30)/10)*100;
            sigma2_n = LPF_cut_margin * (CAL_pow*BLK_pow(2)); 
            mu_H0 = 0;
            N_effect = sig_length/BLK_OS;
            sigma2_H0 = sigma2_n^2/N_effect;
            TH = qfuncinv(0.05)*sqrt(sigma2_H0);
            PFA(pp) = sum(ED_out_H0(:,pp)>TH)/MCtimes;
            PM(pp) = sum(ED_out_H1(:,pp)<TH)/MCtimes;
        end
    else % uses the TH from measureming distribution in H0
        for pp = 1:length(pow_range)
            ED_sort_H = sort(ED_out_H0(:,pp),'descend');
            TH = ED_sort_H(floor(length(ED_sort_H)*0.05));
            PFA(pp) = sum(ED_out_H0(:,pp)>TH)/MCtimes;
            PM(pp) = sum(ED_out_H1(:,pp)<TH)/MCtimes;
        end
    end

    figure
    plot(pow_range,PFA,'o-');hold on
    plot(pow_range,PM,'o-');hold on
    legend('PFA','PM')
end
%% Waveform (for debug)
% BLK BB waveform (original)
% figure;
% plot(linspace(0,200,length(BLK_BB(1:end-2e4,1))),1*IP_mag*BLK_BB(1:end-2e4,1),'linewidth',2);
% hold on;
% % BLK BB waveform (square downconverting and PN despreading; upper arm)
% plot(linspace(0,200,length(LPF_out_H1_u)),LPF_out_H1_u,'--','linewidth',2);
% hold on;
% % BLK BB waveform (square downconverting and PN despreading; lower arm)
% plot(linspace(0,200,length(LPF_out_H1_l)),LPF_out_H1_l,'--','linewidth',2);
% legend('Original','Upper Arm','Lower Arm')

%% PSD plot (for system with freq. offset b/w PN and BLK)
data_num = 1.5e4;
YDATA = fft(DSP_in_H0_u);
figure
subplot(521)
semilogy(linspace(-1,1,length(YDATA)),[abs(YDATA(length(YDATA)/2+1:end));abs(YDATA(1:length(YDATA)/2))])
xlabel('Freq [pi]')
ylabel('PSD')
grid on
title('ADC output (B1 off, B2 on)')
YDATA = fft(DSP_in_H1_u);
subplot(522)
semilogy(linspace(-1,1,length(YDATA)),[abs(YDATA(length(YDATA)/2+1:end));abs(YDATA(1:length(YDATA)/2))])
xlabel('Freq [pi]')
ylabel('PSD')
grid on
title('ADC output (B1 on, B2 on)')


YDATA = fft((DSP_in_H0_u.*carrier_OS));
subplot(523)
semilogy(linspace(-1,1,length(YDATA)),[abs(YDATA(length(YDATA)/2+1:end));abs(YDATA(1:length(YDATA)/2))])
xlabel('Freq [pi]')
ylabel('PSD')
grid on
title('Shifts BLK\timesPN to DC (B1 off, B2 on)')
YDATA = fft((DSP_in_H1_u.*carrier_OS));
subplot(524)
semilogy(linspace(-1,1,length(YDATA)),[abs(YDATA(length(YDATA)/2+1:end));abs(YDATA(1:length(YDATA)/2))])
xlabel('Freq [pi]')
ylabel('PSD')
grid on
title('Shifts BLK\timesPN to DC (B1 on, B2 on)')

YDATA = fft(sig_analogBB_H0_u);
subplot(525)
semilogy(linspace(-1,1,length(YDATA)),[abs(YDATA(length(YDATA)/2+1:end));abs(YDATA(1:length(YDATA)/2))])
xlabel('Freq [pi]')
ylabel('PSD')
grid on
title('Reject Self-Mixing Terms (B1 off, B2 on)')
YDATA = fft(sig_analogBB_H1_u);
subplot(526)
semilogy(linspace(-1,1,length(YDATA)),[abs(YDATA(length(YDATA)/2+1:end));abs(YDATA(1:length(YDATA)/2))])
xlabel('Freq [pi]')
ylabel('PSD')
grid on
title('Reject Self-Mixing Terms (B1 on, B2 on)')

data_num = 990;
YDATA = fft(corr_in_H0_u(1:data_num));
subplot(527)
semilogy(linspace(-1,1,length(YDATA)),[abs(YDATA(length(YDATA)/2+1:end));abs(YDATA(1:length(YDATA)/2))])
xlabel('Freq [pi]')
ylabel('PSD')
grid on
title('Down-Sample to PN Rate (B1 off, B2 on)')
YDATA = fft(corr_in_H1_u(1:data_num));
subplot(528)
semilogy(linspace(-1,1,length(YDATA)),[abs(YDATA(length(YDATA)/2+1:end));abs(YDATA(1:length(YDATA)/2))])
xlabel('Freq [pi]')
ylabel('PSD')
grid on
title('Down-Sample to PN Rate (B1 on, B2 on)')

YDATA = fft(corr_out_H0_u(1:data_num));
subplot(529)
semilogy(linspace(-1,1,length(YDATA)),[abs(YDATA(length(YDATA)/2+1:end));abs(YDATA(1:length(YDATA)/2))])
xlabel('Freq [pi]')
ylabel('PSD')
grid on
title('After PN1 Despreading (B1 off, B2 on)')
ylim([1e-3,1e-1])

YDATA = fft(corr_out_H1_u(1:data_num));
subplot(5,2,10)
semilogy(linspace(-1,1,length(YDATA)),[abs(YDATA(length(YDATA)/2+1:end));abs(YDATA(1:length(YDATA)/2))])
xlabel('Freq [pi]')
ylabel('PSD')
grid on
title('After PN1 Despreading (B1 on, B2 on)')
ylim([1e-3,1e-1])
