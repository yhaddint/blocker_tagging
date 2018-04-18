% The Rapp approximation model should work in baseband. But it is not clear
% whether it can be used to simulate RF-self mixing of square-law devices

% ------ script control parameter --------
clear;clc;
rng(2)
plot_PSD = 1;
plot_rapp = 0;
MCtimes = 1e2;

% ------ system parameters ---------
sig_length = 2e6;
BLK_num = 2;
CAL_pow = 0.002;
BLK_pow = 0.1; % Blocker power as compared to PN power
freq = [0.2,0.25]; % Center frequency of band (GHz)
mag_H1 = [1, 1];
mag_H0 = [0, 1];
fs = 10; % Conceptual Nyquist rate for simulator (GHz)
PN_BW = 0.005; % BW of calibration signal (5MHz)
PN_OS = fs/PN_BW; % Oversampling ratio for calibration signal
BLK_BW = 0.0002; % BW of blocker (200KHz) 
BLK_OS = fs/BLK_BW; % Oversampling ratio for blocker signal

% Lower-sampling rate blocker sample sequences
[ BLK_symb_temp ] = real(get_SC_waveform( 16, sig_length/10, BLK_OS/1e3, BLK_num ));
BLK_symb = BLK_symb_temp./norm(BLK_symb_temp,'fro')*sqrt(sig_length/10*BLK_num);

%% --------- Monte Carlo Sim --------------


for MCindex = 1:MCtimes
    clc
    fprintf('Iteration %d\n',MCindex);
    
    for bb=1:BLK_num
%         BLK_symb(:,bb) = sqrt(1/5) * (randi(4,sig_length/BLK_OS,1)*2-5);
%         BLK_symb(:,bb) = sqrt(1/1) * (randi(2,sig_length/BLK_OS,1)*2-3);
%         BLK_BB(:,bb) = sqrt(0.1)*kron( BLK_symb(:,bb), ones(1e3,1));
        start_index = randi(length(BLK_symb)-sig_length/1e3)-1;
        BLK_capture(:,bb) = BLK_symb(start_index+1:start_index+sig_length/1e3,bb);
        BLK_BB(:,bb) = sqrt(BLK_pow)*interp1(linspace(0,sig_length-1,sig_length/1e3), BLK_capture(:,bb), linspace(0,sig_length-1,sig_length));
        
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
        sig_H0_u(:,bb) = sqrt(CAL_pow) *carrier(:,bb).*sig_BB_H0_u(:,bb);
        sig_H0_l(:,bb) = sqrt(CAL_pow) *carrier(:,bb).*sig_BB_H0_l(:,bb);
        sig_H1_u(:,bb) = sqrt(CAL_pow) *carrier(:,bb).*sig_BB_H1_u(:,bb);
        sig_H1_l(:,bb) = sqrt(CAL_pow) *carrier(:,bb).*sig_BB_H1_l(:,bb);
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
    sig_analogBB_H0_u = (temp_H0_u(2001:end)-CAL_pow/2*(BLK_num))*sqrt(1/CAL_pow);
    sig_analogBB_H0_l = (temp_H0_l(2001:end)-CAL_pow/2*(BLK_num))*sqrt(1/CAL_pow);
    sig_analogBB_H1_u = (temp_H1_u(2001:end)-CAL_pow/2*(BLK_num))*sqrt(1/CAL_pow);
    sig_analogBB_H1_l = (temp_H1_l(2001:end)-CAL_pow/2*(BLK_num))*sqrt(1/CAL_pow);

    % ------ Sampling with aux ADC (~same rate as baseband PN) --------
    DSP_in_H0_u = downsample(sig_analogBB_H0_u(PN_OS/2:end),PN_OS);
    DSP_in_H0_l = downsample(sig_analogBB_H0_l(PN_OS/2:end),PN_OS);
    DSP_in_H1_u = downsample(sig_analogBB_H1_u(PN_OS/2:end),PN_OS);
    DSP_in_H1_l = downsample(sig_analogBB_H1_l(PN_OS/2:end),PN_OS);
    
    % ---- Correlation ----
   
    corr_out_H0_u = DSP_in_H0_u.*PN_u(1:end-1,1);
    corr_out_H0_l = DSP_in_H0_l.*PN_l(1:end-1,1);
    corr_out_H1_u = DSP_in_H1_u.*PN_u(1:end-1,1);
    corr_out_H1_l = DSP_in_H1_l.*PN_l(1:end-1,1);
    
    % ---- LPF ----
    LPF_cut_margin = BLK_BW/PN_BW*1.25;
    bhi_LPF = fir1(40,LPF_cut_margin,'low');
    LPF_out_temp_H0_u = filter(bhi_LPF,1,corr_out_H0_u);
    LPF_out_H0_u = LPF_out_temp_H0_u(21:end);
    LPF_out_temp_H0_l = filter(bhi_LPF,1,corr_out_H0_l);
    LPF_out_H0_l = LPF_out_temp_H0_l(21:end);
    LPF_out_temp_H1_u = filter(bhi_LPF,1,corr_out_H1_u);
    LPF_out_H1_u = LPF_out_temp_H1_u(21:end);
    LPF_out_temp_H1_l = filter(bhi_LPF,1,corr_out_H1_l);
    LPF_out_H1_l = LPF_out_temp_H1_l(21:end);
    
    % ----- Cross-Multiplication and ED -------
    ED_out_H0(MCindex) = mean(LPF_out_H0_u.*LPF_out_H0_l);
    ED_out_H1(MCindex) = mean(LPF_out_H1_u.*LPF_out_H1_l);
end

%% Evaluation
figure(99)
[b,a] = ksdensity(ED_out_H0);
plot(a,b,'linewidth',2);hold on
[b,a] = ksdensity(ED_out_H1);
plot(a,b,'linewidth',2);hold on
grid on
xlabel('ED output');
ylabel('PDF');
legend('H0','H1')
% detection statistic in H0; Based on 2016/09 slides page 20
sigma2_n = LPF_cut_margin * (CAL_pow*BLK_pow); 
mu_H0 = 0;
N_effect = sig_length/BLK_OS;
sigma2_H0 = sigma2_n^2/N_effect;
x = linspace(mu_H0 - 4*sqrt(sigma2_H0), mu_H0 + 4*sqrt(sigma2_H0), 1e3);
y = normpdf(x,mu_H0,sqrt(sigma2_H0));
figure(99)
plot(x,y,'--','linewidth',2);hold on
% detection statistic in H; Based on 2016/09 slides page 20
sigma2_s = CAL_pow*BLK_pow*0.85; 
mu_H1 = sigma2_s;
N_effect = sig_length/BLK_OS;
sigma2_H1 = (2*sigma2_s^2 + 2*sigma2_s*sigma2_n + sigma2_n^2)/N_effect;
% sigma2_H1 = (2*sigma2_s^2)/N_effect;

x = linspace(mu_H1 - 4*sqrt(sigma2_H1), mu_H1 + 4*sqrt(sigma2_H1), 1e3);
y = normpdf(x,mu_H1,sqrt(sigma2_H1));
figure(99)
plot(x,y,'--','linewidth',2);hold on
%%
plot_PSD = 1;
if plot_PSD
    figure
    [ x_freq, PSD_actual ] = get_PSD( corr_out_H1_u, 998, 40 );
    subplot(211)
    plot(x_freq, PSD_actual);hold on
    xlim([-20,20])
    grid on
    xlabel('Freq.')
    ylabel('PSD')
    [ x_freq, PSD_actual ] = get_PSD( corr_out_H1_l, 998, 40 );
    subplot(212)
    plot(x_freq, PSD_actual);hold on
    xlim([-20,20])
    grid on
    xlabel('Freq.')
    ylabel('PSD')

end

