

% ------ script control parameter --------
clear;clc;warning off
rng(2)                              % Random seed
plot_PSD = 1;                       % Whether plot PSD figure
plot_rapp = 0;                      % Whether plot IP/OP relationship in LNA 
MCtimes = 1e2;                    % Number of Monte Carlo trials

% ------ system parameters ---------
T_win = 200e-6;                     % Time window to capture signal in [s], i.e., 200us 
fs = 5;                             % RF Nyquist rate for simulator [GHz]
PN_BW = 0.005;                      % BW of PN signal in [GHz], i.e.,5MHz
PN_OS = fs/PN_BW;                   % Oversampling ratio for PN signal
BLK_BW = 0.0002;                    % BW of blocker in [GHz], i.e., 200KHz
BLK_OS = fs/BLK_BW;                 % Oversampling ratio for blocker signal
sig_length = T_win*fs*1e9;          % Number of sample (conceptual sample in RF domain)
BLK_num = 2;
% CAL_pow = 0.04;                   % Power of calibration (PN) signal (to model saturation)
BLK_pow = [0.1 0.1];                % Blocker power (relative to PN power)
freq = [0.2,0.5];                   % Center frequency of blockers [GHz]
mag_H1 = [1, 1];                    % Enable/disable blockers in hypothesis testing
mag_H0 = [0, 1];                    % Enable/disable blockers in hypothesis testing
resistance = 50;                    % default resistance [om]

% ---- Sample Sequences of Blockers -------
% Sim. uses 2 stages - 1) 50x oversampling with pulse shaping function
%                      2) 500x interpolation without pulse shaping
MQAM = 16;                              % Modulation level, 16QAM
upsam = 50;                             % Oversampling ratio in stage 1
upsam_addition = BLK_OS/upsam;          % Requires interpolation in stage 2
BLK_symb_temp = real( get_SC_waveform( MQAM, 1e5, upsam, BLK_num ) );
BLK_symb = BLK_symb_temp./norm(BLK_symb_temp,'fro')*sqrt(sig_length/10*BLK_num);

% ---- Analog Anti-aliasing Filter -------
fcutoff = 0.015;                          % Fcutoff is around 15MHz
fcutoff_DSP = fcutoff/fs*2;               % Effective "digital" frequency in RF sim.
AA_taps = 2000;                           % Num of Tap (just for RF sim purpose)
bhi_AA = fir1(AA_taps,fcutoff_DSP,'low'); % Fcutoff is around 15MHz

% ---- Filter to reject PN ------
LPF_cut_margin = BLK_BW/PN_BW*1.25;       % LPF in DSP (keep blocker & reject PN spreaded blockers)
bhi_LPF = fir1(40,LPF_cut_margin,'low');  % Filter design

pow_range = [2.5,0,-2.5,-5,-7.5,-10,-15,-25,-35]; % Test various IP power
% pow_range = [-20];                                % Total input power range in [dBm] (BLK and PN)

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
            sig_BB_H0_u(:,bb) = mag_H0(bb) * BLK_BB(:,bb) + CAL_BB_u(:,bb);
            sig_BB_H0_l(:,bb) = mag_H0(bb) * BLK_BB(:,bb) + CAL_BB_l(:,bb);
            sig_BB_H1_u(:,bb) = mag_H1(bb) * BLK_BB(:,bb) + CAL_BB_u(:,bb);
            sig_BB_H1_l(:,bb) = mag_H1(bb) * BLK_BB(:,bb) + CAL_BB_l(:,bb);

            carrier(:,bb) = cos(pi*2*(freq(bb)/fs)*(0:sig_length-1).');

            % ---- RF signals ------
            sig_H0_u(:,bb) = IP_mag * carrier(:,bb) .* sig_BB_H0_u(:,bb);
            sig_H0_l(:,bb) = IP_mag * carrier(:,bb) .* sig_BB_H0_l(:,bb);
            sig_H1_u(:,bb) = IP_mag * carrier(:,bb) .* sig_BB_H1_u(:,bb);
            sig_H1_l(:,bb) = IP_mag * carrier(:,bb) .* sig_BB_H1_l(:,bb);
            
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
        DC_u = (rand*2-1)*1e-4;
        DC_l = (rand*2-1)*1e-4;
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
        
        % ----- Digital compensate signals from PN self-mixing --------
        sig_analogBB_H0_u = (temp_H0_u(AA_shift:end)...
                                -(IP_mag^2)/2*(BLK_num))*sqrt(1/IP_mag^2);
        sig_analogBB_H0_l = (temp_H0_l(AA_shift:end)...
                                -(IP_mag^2)/2*(BLK_num))*sqrt(1/IP_mag^2);
        sig_analogBB_H1_u = (temp_H1_u(AA_shift:end)...
                                -(IP_mag^2)/2*(BLK_num))*sqrt(1/IP_mag^2);
        sig_analogBB_H1_l = (temp_H1_l(AA_shift:end)...
                                -(IP_mag^2)/2*(BLK_num))*sqrt(1/IP_mag^2);

        % ------ Sampling with aux ADC (~same rate as baseband PN) --------
        DSP_in_H0_u = downsample( sig_analogBB_H0_u(PN_OS/2:end), PN_OS );
        DSP_in_H0_l = downsample( sig_analogBB_H0_l(PN_OS/2:end), PN_OS );
        DSP_in_H1_u = downsample( sig_analogBB_H1_u(PN_OS/2:end), PN_OS );
        DSP_in_H1_l = downsample( sig_analogBB_H1_l(PN_OS/2:end), PN_OS );

        % ---- Correlation (Assuming Perfect Timing b/w LNA IP/OP) ----
        corr_out_H0_u = DSP_in_H0_u.*PN_u(1:length(DSP_in_H0_u),1);
        corr_out_H0_l = DSP_in_H0_l.*PN_l(1:length(DSP_in_H0_l),1);
        corr_out_H1_u = DSP_in_H1_u.*PN_u(1:length(DSP_in_H1_u),1);
        corr_out_H1_l = DSP_in_H1_l.*PN_l(1:length(DSP_in_H1_l),1);

        % ---- LPF to Reject Spreaded Signals ----
        LPF_out_temp_H0_u = filter(bhi_LPF,1,corr_out_H0_u);
        LPF_out_H0_u = LPF_out_temp_H0_u(21:end);
        LPF_out_temp_H0_l = filter(bhi_LPF,1,corr_out_H0_l);
        LPF_out_H0_l = LPF_out_temp_H0_l(21:end);
        LPF_out_temp_H1_u = filter(bhi_LPF,1,corr_out_H1_u);
        LPF_out_H1_u = LPF_out_temp_H1_u(21:end);
        LPF_out_temp_H1_l = filter(bhi_LPF,1,corr_out_H1_l);
        LPF_out_H1_l = LPF_out_temp_H1_l(21:end);

        % ----- Cross-Multiplication and ED -------
        ED_out_H0(MCindex,pp) = mean(LPF_out_H0_u.*LPF_out_H0_l);
        ED_out_H1(MCindex,pp) = mean(LPF_out_H1_u.*LPF_out_H1_l);
        
    end % end of Monte Carlo loop
end % end of IP power sweeping

%% Evaluation the detection statistic (for debug)
pp=1;
figure(99)
[b,a] = ksdensity(ED_out_H0(:,pp));
plot(a,b,'linewidth',2);hold on
[b,a] = ksdensity(ED_out_H1(:,pp));
plot(a,b,'linewidth',2);hold on
grid on
xlabel('ED output');
ylabel('PDF');
legend('H0','H1')
% detection statistic in H0; Based on 2016/09 slides page 20
CAL_pow(pp) = 10^((pow_range(pp)-30)/10)*100;
sigma2_n = LPF_cut_margin * (CAL_pow(pp)*BLK_pow(2)); 
mu_H0 = 0;
N_effect = sig_length/BLK_OS;
sigma2_H0 = sigma2_n^2/N_effect;
x = linspace(mu_H0 - 4*sqrt(sigma2_H0), mu_H0 + 4*sqrt(sigma2_H0), 1e3);
y = normpdf(x,mu_H0,sqrt(sigma2_H0));
figure(99)
plot(x,y,'--','linewidth',2);hold on
% detection statistic in H; Based on 2016/09 slides page 20
sigma2_s = CAL_pow(pp)*BLK_pow(1)*0.85; 
mu_H1 = sigma2_s;
N_effect = sig_length/BLK_OS;
sigma2_H1 = (2*sigma2_s^2 + 2*sigma2_s*sigma2_n + sigma2_n^2)/N_effect;
% sigma2_H1 = (2*sigma2_s^2)/N_effect;

x = linspace(mu_H1 - 4*sqrt(sigma2_H1), mu_H1 + 4*sqrt(sigma2_H1), 1e3);
y = normpdf(x,mu_H1,sqrt(sigma2_H1));
figure(99)
plot(x,y,'--','linewidth',2);hold on
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
for pp=1:length(pow_range)
    CAL_pow = 10^((pow_range(pp)-30)/10)*100;
    sigma2_n = LPF_cut_margin * (CAL_pow*BLK_pow(2)); 
    mu_H0 = 0;
    N_effect = sig_length/BLK_OS;
    sigma2_H0 = sigma2_n^2/N_effect;
    TH = qfuncinv(0.05)*sqrt(sigma2_H0);
    PFA(pp) = sum(ED_out_H0(:,pp)>TH)/MCtimes;
    PM(pp) = sum(ED_out_H1(:,pp)<TH)/MCtimes;
end
figure
plot(pow_range,PFA,'o-');hold on
plot(pow_range,PM,'o-');hold on
legend('PFA','PM')
%% Waveform (for debug)
% BLK BB waveform (original)
figure;
plot(linspace(0,200,length(BLK_BB(1:end-2e4,1))),sqrt(CAL_pow)*BLK_BB(1:end-2e4,1),'linewidth',2);
hold on;
% BLK BB waveform (square downconverting and PN despreading; upper arm)
plot(linspace(0,200,length(LPF_out_H1_u)),LPF_out_H1_u,'--','linewidth',2);
hold on;
% BLK BB waveform (square downconverting and PN despreading; lower arm)
plot(linspace(0,200,length(LPF_out_H1_l)),LPF_out_H1_l,'--','linewidth',2);