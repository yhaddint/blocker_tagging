% ------ script control parameter --------
clear;clc;warning off
rng(2)                              % Random seed

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
IPpow = -10;                        % Total input power in [dBm] (BLK and PN)
IP_mag = dBm2V( IPpow, resistance ); % Input magnitude in [Volt]
freq_offset = 10e6;                         % offset b/w BLK and PN

% ---- Sample Sequences of Blockers -------
% Sim. uses 2 stages - 1) 50x oversampling with pulse shaping function
%                      2) 500x interpolation without pulse shaping
MQAM = 16;                              % Modulation level, 16QAM
upsam = 50;                             % Oversampling ratio in stage 1
upsam_addition = BLK_OS/upsam;          % Requires interpolation in stage 2
BLK_symb_temp = real( get_SC_waveform( MQAM, 1e5, upsam, BLK_num ) );
BLK_symb = BLK_symb_temp./norm(BLK_symb_temp,'fro')*sqrt(sig_length/10*BLK_num);

  
% Waveform per band
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
    sig_BB_H0_u(:,bb) = mag_H0(bb) * BLK_BB(:,bb);
    sig_BB_H0_l(:,bb) = mag_H0(bb) * BLK_BB(:,bb);
    sig_BB_H1_u(:,bb) = mag_H1(bb) * BLK_BB(:,bb);
    sig_BB_H1_l(:,bb) = mag_H1(bb) * BLK_BB(:,bb);

    carrier_BLK(:,bb) = cos(pi*2*(freq(bb)/fs)*(0:sig_length-1).');
    carrier_PN(:,bb) = cos(pi*2*((freq(bb) + freq_offset)/fs)*(0:sig_length-1).');


    % ---- RF signals ------
    sig_H0_u(:,bb) = IP_mag * (carrier_BLK(:,bb) .* sig_BB_H0_u(:,bb)...
                              +carrier_PN(:,bb) .* CAL_BB_u(:,bb) );
                          
    sig_H0_l(:,bb) = IP_mag * (carrier_BLK(:,bb) .* sig_BB_H0_l(:,bb)...
                              +carrier_PN(:,bb) .* CAL_BB_l(:,bb) );
                          
    sig_H1_u(:,bb) = IP_mag * (carrier_BLK(:,bb) .* sig_BB_H1_u(:,bb)...
                              +carrier_PN(:,bb) .* CAL_BB_u(:,bb) );
                          
    sig_H1_l(:,bb) = IP_mag * (carrier_BLK(:,bb) .* sig_BB_H1_l(:,bb)...
                              +carrier_PN(:,bb) .* CAL_BB_l(:,bb) );

end

% ---- Input RF Waveform of LNA -------
sig_in_H0_u = sum(sig_H0_u,2);
sig_in_H0_l = sum(sig_H0_l,2);
sig_in_H1_u = sum(sig_H1_u,2);
sig_in_H1_l = sum(sig_H1_l,2);
RF_waveform = [sig_in_H0_u,...
               sig_in_H0_l,...
               sig_in_H1_u,...
               sig_in_H1_l];
           
IP_RMS_pow_H0 = V2dBm( sqrt(mean(abs(sig_in_H0_u).^2)) );
IP_RMS_pow_H1 = V2dBm( sqrt(mean(abs(sig_in_H1_u).^2)) );

% --- print iteration number ----
fprintf('IP power (requested) %4.2f dBm\n', IPpow);
fprintf('IP power (actual in H0) %4.2f dBm\n', IP_RMS_pow_H0);
fprintf('IP power (actual in H1) %4.2f dBm\n', IP_RMS_pow_H1);

% save data
filename = 'RF_waveform.mat';
save(filename, 'RF_waveform');



