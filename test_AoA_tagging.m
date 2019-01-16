clear;clc;
%  11/13/2018
%  path tagging idea. Concept adapts from single arm blocker tagging

clear;clc; warning off %clf;close all


% rand('seed',3)
% ---------- system parameters ----------
path_num = 2;                       % Num. of AoA sources
M = 128;                            % Num. of symbol (PSS/Golay) in each sounding BF
Ns = 40;                            % Num. of sample
SampleNumberAve = Ns;               % Oversampling ratio in Tx PN v.s. signal
stat_num = 2e2;                     % Control length of received signal
P = M*Ns;                           %
L = M*SampleNumberAve*stat_num;     % control singal length
MCtimes = 1e2;                      % Monte Carlo simulation number
Ntr = 16;                           % TRx array size for AoA search
NB = Ntr*2;                           % Num. of inserted signals
SNR_dB = 0;
SNR = 10^(SNR_dB/10);
% ---------- Signal Waveform Generation ----------
for ii = 1:path_num
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
    upsam = SampleNumberAve;
    symbols = fix(L*2/upsam);
    
    clearvars data temp_data
    hmod = modem.qammod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
%     unnormalized = real(temp_data(end-L+1:end));
%     sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     unnormalized = imag(temp_data(end-L+1:end));
%     sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     sig_i(:,ii) = zeros(size(sig_r(:,ii)));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
    sig(:,ii) = temp_data(end-L+1:end)./sqrt(temp_data(end-L+1:end)'*temp_data(end-L+1:end)/L);
end

%     ii = 1;
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
%     upsam = 128;
%     L2 = M*upsam*stat_num;
%     symbols = fix(L2*2/upsam);
%     
%     clearvars data temp_data
%     hmod = modem.qammod('M', 4, 'InputType', 'integer');
%     hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
%     hpulse = design(hdesign);
%     data = randi(4,symbols,1)-1;
%     data = modulate(hmod, data);
%     data = upsample(data,upsam);
%     temp_data = conv(data,hpulse.Numerator);
%     unnormalized = real(temp_data(end-L+1:end));
%     sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     unnormalized = imag(temp_data(end-L+1:end));
%     sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     sig_i(:,ii) = zeros(size(sig_r(:,ii)));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
%     sig_other = temp_data(end-L+1:end)./sqrt(temp_data(end-L+1:end)'*temp_data(end-L+1:end)/L);


%%
other_pow = 1e-4;
AoA_pow0 = ones(path_num,1)*other_pow;
AoA_pow1 = ones(path_num,1)*other_pow;


AoA_pow1(1:2) = [1,1e-1]';
AoA_pow0(1:2) = [1e-3,1e-3]';

AoA_angle = [-24,30].'/180*pi; % AoA in propagation direction
AoA_lambda = pi*sin(AoA_angle); % AoA in the phase of wavelength

psi_angle = linspace(-90,90,NB).'/180*pi;
psi_lambda = pi*sin(psi_angle);
%% Monte Carlo Simulation

for MCindex=1:MCtimes
    clc; fprintf('Iteration %d out of %d\n', MCindex, MCtimes)

    % ------ Generate PN signal ---------
    CAL = (randi(2,P,NB)*2-3)/sqrt(NB);

    % ------ Generate signal ---------
    start_point = randi(L-P-1,path_num,1);
    for ii=1:path_num
        sig_cache(:,ii) = sig(start_point(ii,1):start_point(ii,1)+P-1,randi(path_num));
    end
%     sig_cache(:,4) = sig_other(start_point(ii,1):start_point(ii,1)+P-1,1);
    
    % ------ Zero Initialization of vec/mtx ---------
    sig_mix0_part1 = zeros(P,Ntr);
    sig_mix1_part1 = zeros(P,Ntr);
%     sig_mix0_part2 = zeros(P,Ntr);
%     sig_mix1_part2 = zeros(P,Ntr);   
%     sig_mix0_part3 = zeros(P,Ntr);
%     sig_mix1_part3 = zeros(P,Ntr);
    sig_mix0 = zeros(P,Ntr);
    sig_mix1 = zeros(P,Ntr);
    
    % Unknown channel gain in each path (equal power for now)
    PhaseShift = exp(1j*rand(1,path_num));
    
    % -------- Baseband Simulation of RF Nonlinear Terms (Sig x PN)-------
    for nn=1:Ntr
        for ii=1:path_num
            for bb=1:NB
    %             carrier = exp(1j*((ii-jj)/32*2*pi*(1:P)'))*PhaseShift(ii);
                sig_mix0_part1(:,nn) = sig_mix0_part1(:,nn) +...
                    sqrt(AoA_pow0(ii))*CAL(:,bb).* ...
                    real(sig_cache(:,ii)*exp(1j*(nn-1)*(AoA_lambda(ii)-psi_lambda(bb))));
                sig_mix1_part1(:,nn) = sig_mix1_part1(:,nn) +...
                    sqrt(AoA_pow1(ii))*CAL(:,bb).* ...
                    real(sig_cache(:,ii)*exp(1j*(nn-1)*(AoA_lambda(ii)-psi_lambda(bb))));
            end
        end
    end
    
    % ----------- Noise Generation -----------
    noise_mag = sqrt(norm(sig_mix1_part1,'fro')^2/Ntr/P/SNR);
    AWGN = (randn(P,Ntr) + 1j*randn(P,Ntr))/sqrt(2)*noise_mag;
    
    % -------- Baseband Simulation of RF Nonlinear Terms (All)-------
    sig_mix0 = sum(sig_mix0_part1 + AWGN,2); % Spatial filtering (sum for now)
    sig_mix1 = sum(sig_mix1_part1 + AWGN,2); % Spatial filtering (sum for now)
    
    % --------- DSP Tagging Algorithm -----------
    bhi = fir1(200,1/Ns,'low');
    for bb=1:NB
        sig_result1(:,bb) = sqrt(NB)*filter(bhi,1,sig_mix1(1:P).*CAL(1:P,bb));
        ED_results1(bb,MCindex) = norm(sig_result1(101:end,bb));
        
        sig_result0(:,bb) = sqrt(NB)*filter(bhi,1,sig_mix0(1:P).*CAL(1:P,bb));
        ED_results0(bb,MCindex) = norm(sig_result0(101:end,bb));
    end
end

%% FFT test
% bb_true =7;
% xdata = linspace(-pi,pi,P);
% ytemp = abs(fft(sig_mix1(1:P).*CAL(1:P,bb_true)));
% ytemp2 = abs(fft(sig_result1(:,bb_true)));
% 
% figure
% plot(xdata,20*log10([ytemp(P/2:end);ytemp(1:P/2-1)]));hold on
% plot(xdata,20*log10([ytemp2(P/2:end);ytemp2(1:P/2-1)]));hold on
% xlabel('Freq. [rad]')
% grid on
%% plot
figure;
plot(psi_angle/pi*180, 20*log10(mean(ED_results1,2)),'-o');hold on

for bb=1:NB
    temp = sort(20*log10(ED_results1(bb,:)),'descend');
    upper(bb) = temp(floor(MCtimes*0.1));
    lower(bb) = temp(floor(MCtimes*0.9));
    plot(ones(2,1)*psi_angle(bb)/pi*180, [lower(bb) upper(bb)] ,'k-');hold on
end


% plot(psi_angle/pi*180, upper,'.');hold on
% plot(psi_angle/pi*180, lower,'.');hold on

% plot(psi_angle/pi*180, 10*log10(mean(ED_results0,2)),'--x');hold on
xlabel('Angle Candidates [deg]')
ylabel('Score [dB]')
grid on
dim = [.6 .6 .3 .3];
str1 = ['True AoA1 = ' num2str(AoA_angle(1)/pi*180) ' deg'];
str2 = ['True AoA2 = ' num2str(AoA_angle(2)/pi*180) ' deg'];
str = {str1,str2};
annotation('textbox',dim,'String',str,'FitBoxToText','on');


