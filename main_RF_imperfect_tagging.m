% The Rapp approximation model should work in baseband. But it is not clear
% whether it can be used to simulate RF-self mixing of square-law devices

% ------ script control parameter --------
clear;clc;
rng(2)
plot_PSD = 1;
plot_rapp = 0;
MCtimes = 1e2;

% ------ system parameters ---------
sig_length = 1e6;
BLK_num = 2;
CAL_pow = 0.002;
freq = [0.2,0.25]; % Center frequency of band (GHz)
mag_H1 = [1, 1];
mag_H0 = [0, 1];
fs = 10; % Conceptual Nyquist rate for simulator (GHz)
PN_BW = 0.005; % BW of calibration signal (5MHz)
PN_OS = fs/PN_BW; % Oversampling ratio for calibration signal
BLK_BW = 0.0002; % BW of blocker (200KHz) 
BLK_OS = fs/BLK_BW; % Oversampling ratio for blocker signal

%% --------- Monte Carlo Sim --------------
for MCindex = 1:MCtimes
    clc
    fprintf('Iteration %d\n',MCindex);
    for bb=1:BLK_num
        BLK_BB(:,bb) = 0.05*kron(randi(4,sig_length/BLK_OS,1)*2-5,ones(BLK_OS,1));
        % ---- PN sequence -------
        PN_u(:,bb) = randi(2,sig_length/PN_OS,1)*2-3;
        PN_l(:,bb) = randi(2,sig_length/PN_OS,1)*2-3;
        
        % ---- PN sequence -------
        CAL_BB_u(:,bb) = kron(PN_u(:,bb),ones(PN_OS,1));
        CAL_BB_l(:,bb) = kron(PN_l(:,bb),ones(PN_OS,1));
        
        % ---- In-band signal (baseband version) -------
        sig_BB_H0_u(:,bb) = mag_H0(bb) * BLK_BB(:,bb) + CAL_BB_u(:,bb);
        sig_BB_H0_l(:,bb) = mag_H0(bb) * BLK_BB(:,bb) + CAL_BB_l(:,bb);
        sig_BB_H1_u(:,bb) = mag_H1(bb) * BLK_BB(:,bb) + CAL_BB_u(:,bb);
        sig_BB_H1_l(:,bb) = mag_H1(bb) * BLK_BB(:,bb) + CAL_BB_l(:,bb);
        
        carrier(:,bb) = cos(pi*2*(freq(bb)/fs)*(0:sig_length-1).');
        
        % ---- RF signals ------
        sig_H0_u(:,bb) = sqrt(CAL_pow*2) *carrier(:,bb).*sig_BB_H0_u(:,bb);
        sig_H0_l(:,bb) = sqrt(CAL_pow*2) *carrier(:,bb).*sig_BB_H0_l(:,bb);
        sig_H1_u(:,bb) = sqrt(CAL_pow*2) *carrier(:,bb).*sig_BB_H1_u(:,bb);
        sig_H1_l(:,bb) = sqrt(CAL_pow*2) *carrier(:,bb).*sig_BB_H1_l(:,bb);
    end
    
    % ----- Rapp model for square-law devices -----
    Asat = 1;
    P_rapp = 1;
    sig_out_H0_u = get_rapp_square(sum(sig_H0_u,2),Asat,P_rapp);
    sig_out_H0_l = get_rapp_square(sum(sig_H0_l,2),Asat,P_rapp);
    sig_out_H1_u = get_rapp_square(sum(sig_H1_u,2),Asat,P_rapp);
    sig_out_H1_l = get_rapp_square(sum(sig_H1_l,2),Asat,P_rapp);
    
    % ---- Anti-Aliasing LPF before aux ADC -----
    fcutoff = 0.015; % Fcutoff is around 15MHz
    fcutoff_DSP = fcutoff/fs*2;
    bhi = fir1(4000,fcutoff_DSP,'low'); % Fcutoff is around 15MHz
    temp_H0_u = filter(bhi,1,sig_out_H0_u);
    temp_H0_l = filter(bhi,1,sig_out_H0_l);
    temp_H1_u = filter(bhi,1,sig_out_H1_u);
    temp_H1_l = filter(bhi,1,sig_out_H1_l);
    sig_analogBB_H0_u = temp_H0_u(2001:end)-CAL_pow*(BLK_num);
    sig_analogBB_H0_l = temp_H0_l(2001:end)-CAL_pow*(BLK_num);
    sig_analogBB_H1_u = temp_H1_u(2001:end)-CAL_pow*(BLK_num);
    sig_analogBB_H1_l = temp_H1_l(2001:end)-CAL_pow*(BLK_num);

    % ------ Sampling with aux ADC (~10MHz) --------
    DSP_in_H0_u = downsample(sig_analogBB_H0_u(1e3:end),2e3);
    DSP_in_H0_l = downsample(sig_analogBB_H0_l(1e3:end),2e3);
    DSP_in_H1_u = downsample(sig_analogBB_H1_u(1e3:end),2e3);
    DSP_in_H1_l = downsample(sig_analogBB_H1_l(1e3:end),2e3);
    
    % Correlation & LPF
    bhi_ED = fir1(40,0.05,'low');
    corr_out_H0 = (DSP_in_H0_u.*PN_u(1:end-1,1)).*(DSP_in_H0_l.*PN_l(1:end-1,1));
    corr_out_H1 = (DSP_in_H1_u.*PN_u(1:end-1,1)).*(DSP_in_H1_l.*PN_l(1:end-1,1));
    ED_out_temp_H0 = filter(bhi_ED,1,corr_out_H0);
    ED_out_H0(MCindex) = sum(abs(ED_out_temp_H0(21:end)).^2);
    ED_out_temp_H1 = filter(bhi_ED,1,corr_out_H1);
    ED_out_H1(MCindex) = sum(abs(ED_out_temp_H1(21:end)).^2);
end

%% Evaluation
figure
[b,a] = ksdensity(ED_out_H0);
plot(a,b);hold on
[b,a] = ksdensity(ED_out_H1);
plot(a,b);hold on
grid on
xlabel('ED output');
ylabel('PDF');
legend('H0','H1')
%%

% figure
% [ x_freq, PSD_actual ] = get_PSD( corr_out_H0, 498, pi );
% subplot(211)
% plot(x_freq, PSD_actual);hold on
% % xlim([0,1e3])
% grid on
% xlabel('Freq')
% ylabel('PSD')
% [ x_freq, PSD_actual ] = get_PSD( corr_out_H1, 498, pi );
% subplot(212)
% plot(x_freq, PSD_actual);hold on
% % xlim([0,1e3])
% grid on
% xlabel('Freq')
% ylabel('PSD')


