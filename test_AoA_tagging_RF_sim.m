clear;clc;
%  11/13/2018
%  path tagging idea. Concept adapts from single arm blocker tagging

clear;clc; warning off %clf;close all


rand('seed',3)
% ---------- system parameters ----------
T_win = 200e-6;                     % Capture signal in a 200us window
fs = 5;                             % Conceptual Nyquist rate for simulator [GHz]
PN_BW = 0.005;                      % BW of calibration signal 5[MHz]
PN_OS = fs/PN_BW;                   % Oversampling ratio for calibration signal
SIG_BW = 0.0002;                    % BW of signal source 200[KHz]
SIG_OS = fs/SIG_BW;                 % Oversampling ratio for signal
sig_length = T_win*fs*1e9;          % Number of sample (conceptual sample in RF domain)

path_num = 2;                       % Num. of AoA sources
M = 128;                            % Num. of symbol (PSS/Golay) in each sounding BF
SampleNumberAve = 40;               % Oversampling ratio in Tx PN v.s. signal
P = 1e4;                            % control singal length
L = 2e4;                            % control singal length
MCtimes = 2e1;                      % Monte Carlo simulation number
Ntr = 16;                           % TRx array size for AoA search
NB = Ntr * 2;                       % Num. of inserted signals
SNR_dB = 0;
SNR = 10^(SNR_dB/10);
freq = 0.2;                         % Freq. for RF sim purpose [unitless]
DAC_range = 3:2:15;
DAC_num = length(DAC_range);

% ---- Lower-sampling rate blocker sample sequences -------
MQAM = 16; % Modulation level
upsam = 50;
upsam_addition = SIG_OS/upsam;
[ SIG_symb_temp ] = real(get_SC_waveform( MQAM, 1e5, upsam, path_num ));
SIG_symb = SIG_symb_temp./norm(SIG_symb_temp,'fro')*sqrt(sig_length/10*path_num);

% ---- Anti-aliasing filter  -------
fcutoff = 0.015; % Fcutoff is around 15MHz
fcutoff_DSP = fcutoff/fs*2;
AA_taps = 2000;
bhi_AA = fir1(AA_taps,fcutoff_DSP,'low'); % Fcutoff is around 15MHz

% ---- Filter to reject PN ------
LPF_cut_margin = SIG_BW/PN_BW*1.25;
bhi_LPF = fir1(40,LPF_cut_margin,'low');

bb_true = 8;
%%
AoA_pow(1:2) = [1e-1,0]';
CAL_pow = 10^((-40-30)/10)*100;

AoA_angle = [-33,30].'/180*pi; % AoA in propagation direction
AoA_lambda = pi*sin(AoA_angle); % AoA in the phase of wavelength

psi_angle = linspace(-60,60,NB).'/180*pi;
psi_lambda = pi*sin(psi_angle);
%% Monte Carlo Simulation

for MCidx = 1 : MCtimes
    clc; fprintf('Iteration %d out of %d\n', MCidx, MCtimes)
    
    
    % ---- PN Sequence -------
    for PN_idx = 1 : NB
        PN(:,PN_idx) = randi(2,sig_length/PN_OS,1)*2-3;
        CAL_BB(:,PN_idx) = kron(PN(:,PN_idx),ones(PN_OS,1));
    end
    
    % ------ Generate Signal ---------
    for path_idx = 1 : path_num
        start_index = randi(size(SIG_symb, 1) - sig_length/upsam_addition) - 1;
        SIG_capture(:,path_idx) = SIG_symb(start_index+1:start_index+sig_length/upsam_addition,path_idx);
        time_orig = linspace(0,sig_length-1,sig_length/upsam_addition);
        time_intep = linspace(0,sig_length-1,sig_length);
        SIG_BB(:,path_idx) = sqrt(AoA_pow(path_idx))*interp1(time_orig, SIG_capture(:,path_idx), time_intep);
    end
    
    % -------- Baseband Sim. of RF Nonlinear Terms (Sig + many PN + AWGN)-------
    carrier = exp(1j*pi*2*(freq/fs)*(0:sig_length-1).');
    sig_mix_PN = zeros(sig_length, Ntr);
    PN_per_ant = zeros(sig_length, Ntr);
    for nn = 1:Ntr
        
        % ---- PN with phase shifter, "spatially" inserted -------
        for PN_idx = 1 : NB
            if mod(PN_idx,2)==0
                PN_per_ant(:,nn) = PN_per_ant(:,nn) ...
                    + CAL_BB(:,PN_idx).* exp(1j*(nn-1)*(psi_lambda(PN_idx)));
            end
        end

        % ------ Original received signal -----------
        for path_idx = 1 : path_num
            sig_mix_PN(:,nn) = sig_mix_PN(:,nn) + sqrt(AoA_pow(path_idx))...
                *ones(sig_length,1) .* real(carrier*exp(1j*(nn-1)*(AoA_lambda(path_idx))));
        end
        
    end
        
    for DAC_idx = 1:DAC_num
        DAC_bits = DAC_range(DAC_idx);
       
        for nn = 1:Ntr    
            
            % -------- Apply quantization ---------
            max_mag = max(abs(PN_per_ant(:,nn)));
            PN_quan(:,nn) = DAC_quan( PN_per_ant(:,nn), DAC_bits, max_mag );
            PN_quan_RF(:,nn) = real(PN_quan(:,nn).*carrier);

            % ----------- Noise Generation -----------
    %         noise_mag = sqrt(norm(sig_mix_PN,'fro')^2/Ntr/P/SNR);
    %         AWGN = (randn(P,Ntr) + 1j*randn(P,Ntr))/sqrt(2)*noise_mag;

            % Adding all signal together at each antenna
            sig_RF(:,nn) = sig_mix_PN(:,nn) + PN_quan_RF(:,nn);% + AWGN;
            null_mix_PN(:,nn) = PN_quan_RF(:,nn);

            % ------ Rapp model for square-law devices --------
            Asat = 1;
            P_rapp = 1;
            sig_mix_PN_out(:,nn) = get_rapp_square( sig_RF(:,nn), Asat, P_rapp);
            null_mix_PN_out(:,nn) = get_rapp_square( null_mix_PN(:,nn), Asat, P_rapp);


            % ------ Downsample for BB Processing --------
            AA_shift = AA_taps/2+1; % Delay after AA filter
            temp_out = filter(bhi_AA,1,sig_mix_PN_out(:,nn));
            null_temp_out = filter(bhi_AA,1,null_mix_PN_out(:,nn));

            % ------ Subtract DC Components ----------
            sig_analogBB(:,nn) = temp_out(AA_shift:end);%-mean(temp_out(AA_shift:end));
            null_analogBB(:,nn) = null_temp_out(AA_shift:end);%-mean(temp_out(AA_shift:end));


            % ------ Sampling with aux ADC (~same rate as baseband PN) --------
            RF_PS_in(:,nn) = downsample(sig_analogBB(PN_OS/2:end,nn),PN_OS);
            null_PS_in(:,nn) = downsample(null_analogBB(PN_OS/2:end,nn),PN_OS);
        end

        % ------- RF sum (spatial filtering that picks 0 deg) --------
        DSP_in = sum(RF_PS_in,2);
        null_DSP_in = sum(null_PS_in,2);

        % --------- DSP Tagging Algorithm -----------
        for PN_idx = 1 : NB

            % Correlation with PN inserted for each direction
            corr_out = DSP_in.*PN(1:length(DSP_in),PN_idx);
            null_corr_out = null_DSP_in.*PN(1:length(null_DSP_in),PN_idx);

            % Spatial filtering that keep "tagged" path
            LPF_out_temp = filter(bhi_LPF,1,(corr_out-null_corr_out));
            LPF_out(:,PN_idx) = LPF_out_temp(21:end);

            % Energy detection
            ED_results1(PN_idx, MCidx) = mean(LPF_out(:,PN_idx).*LPF_out(:,PN_idx));

         end

        Peak_ratio(DAC_idx,MCidx) = ED_results1(bb_true, MCidx)/max(ED_results1(setdiff(1:NB,bb_true), MCidx));
        
    end
end

%% peak ratio evaluation
plot(DAC_range, 10*log10(mean(Peak_ratio,2)));
hold on
grid on
xlabel('DAC Quantization [bit]')
ylabel('Peak Ratio [dB]')
%%


figure
plot(psi_angle/pi*180, 10*log10(ED_results1(:,1)))
grid on
%% FFT test

xdata = linspace(-pi,pi,999);
ytemp = abs(fft(DSP_in.*CAL_BB(1:length(DSP_in),bb_true)));
% ytemp2 = abs(fft(sig_result1(:,bb_true)));

figure
plot(xdata,20*log10([ytemp(1e3/2:end);ytemp(1:1e3/2-1)]));hold on
% plot(xdata,20*log10([ytemp2(P/2:end);ytemp2(1:P/2-1)]));hold on
xlabel('Freq. [rad]')
grid on
%% plot
figure;
plot(psi_angle/pi*180, 20*log10(mean(ED_results1,2)),'-o');hold on

% for bb=1:NB
%     temp = sort(20*log10(ED_results1(bb,:)),'descend');
%     upper(bb) = temp(floor(MCtimes*0.1));
%     lower(bb) = temp(floor(MCtimes*0.9));
%     plot(ones(2,1)*psi_angle(bb)/pi*180, [lower(bb) upper(bb)] ,'k-');hold on
% end


% plot(psi_angle/pi*180, upper,'.');hold on
% plot(psi_angle/pi*180, lower,'.');hold on

% plot(psi_angle/pi*180, 10*log10(mean(ED_results0,2)),'--x');hold on
xlabel('Angle Candidates [deg]')
ylabel('Score [dB]')
grid on
% dim = [.6 .6 .3 .3];
% str1 = ['True AoA1 = ' num2str(AoA_angle(1)/pi*180) ' deg'];
% str2 = ['True AoA2 = ' num2str(AoA_angle(2)/pi*180) ' deg'];
% str = {str1,str2};
% annotation('textbox',dim,'String',str,'FitBoxToText','on');


