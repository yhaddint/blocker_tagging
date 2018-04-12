%  01/17/2017
%  There are 7 cases. In case 1-6, i^th blocker is turned off and system aims to tag all
%  blockers. therefore in case i, tagging rate of Ji is false alarm rate
%  and the rest is pd. in case 7, all blockers are turned on.

clear;clc;%clf;close all
rng('default')

blk_num = 6;

% oversampling ratio, which is PN rate v.s. blocker rate
OverSampling = 2;

total_time = 1e-3; % maximum running time is 1ms
Lmax = fix(total_time/(1/6e6/OverSampling));

% Length of PN per period
M = 1600;
N = fix(Lmax/M);
P = M*N;
L = 30*P;

Nrange = fix(P*(0.1:0.1:0.9));

Pc = 10; % power of calibration signal (mW)

% generating complex baseband signal, which is the sig_bb used later
for ii = 1:blk_num
    upsam = OverSampling;
    bhi = fir1(400,1/upsam*0.7,'low');
    symbols = fix(L*2/upsam);
    clearvars data temp_data
    hmod = modem.qammod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
    temp1_data = filter(bhi,1,temp_data);
    sig(:,ii) = temp1_data(end-L+1:end)./sqrt(temp1_data(end-L+1:end)'*temp1_data(end-L+1:end)/L);
end


%%
ReusePN = 1;
ofNum = blk_num+1;
MCtimes = 2e2;
results = zeros(ofNum,blk_num,MCtimes,length(Nrange));
% results_H1 = zeros(ofNum,blk_num,runtimes,length(Nrange));
% blocker power settings
sig_pow_dB = [-3,-6,-9,-12,-15,-18]';

sig_pow_all = kron(ones(1,blk_num+1),10.^(sig_pow_dB/10));
for kk=1:blk_num
    sig_pow_all(kk,kk) = 0;
end

for onoffindex = 1:blk_num+1

    XXX = ['simulating, ',num2str(onoffindex/(blk_num+1),2),' finished'];
    disp(XXX)
    sig_pow = sig_pow_all(:,onoffindex);
    
    for MCindex=1:MCtimes
        
        % generating PN sequence
        if ReusePN
            CALseed = randi(2,M+1,blk_num)*2-3;
            CAL_u0 = kron(ones(N,1),CALseed);
            CAL_u = CAL_u0(1:P,:);
            CALseed = randi(2,M+2,blk_num)*2-3;
            CAL_l0 = kron(ones(N,1),CALseed);
            CAL_l = CAL_l0(1:P,:);
        else
            CAL_u = randi(2,P,blk_num)*2-3;
            CAL_l = randi(2,P,blk_num)*2-3;
        end

        % pilots contains two co-channel PN sequences passing through diode
        d = randi(5e3);
        pilot_u = 0.5*Pc*[zeros(d,1);CALseed(:,1)].*conj([zeros(d,1);CALseed(:,1)])...
            +0.5*Pc*[zeros(d,1);CALseed(:,2)].*conj([zeros(d,1);CALseed(:,2)])...
            +Pc*real([zeros(d,1);CALseed(:,1)].*conj([zeros(d,1);CALseed(:,2)]));
        [peak,dhat] = max(abs(conv(pilot_u,flipud(CALseed(:,1).*CALseed(:,2)))));
        % estimated delay offset that will be used in tagging
        dOS = dhat-size(CALseed(:,1),1);
        
        % Simulating Diode after blockers are allowed to bt input
        sig_mix_u = zeros(2*P,1);
        sig_mix_l = zeros(2*P,1);
        start_point = randi(L-P-1,blk_num,1);
        for ii=1:blk_num
            sig_bb(:,ii) = sig(start_point(ii,1):start_point(ii,1)+P-1,randi(blk_num));
        end
        Phi = exp(1j*rand(1,blk_num));
        for kk=1:blk_num
            for ll=1:blk_num
                if kk==ll
                    sig_mix_u(d+1:d+P) = sig_mix_u(d+1:d+P)...
                        +sqrt(Pc*sig_pow(kk))*real(sig_bb(1:P,kk)*Phi(kk)).*CAL_u(1:P,kk)...
                        +0.5*sig_pow(kk)*(sig_bb(1:P,kk).*conj(sig_bb(1:P,kk)));
                    sig_mix_l(d+1:d+P) = sig_mix_l(d+1:d+P)...
                        +sqrt(Pc*sig_pow(kk))*real(sig_bb(1:P,kk)*Phi(kk)).*CAL_l(1:P,kk)...
                        +0.5*sig_pow(kk)*(sig_bb(1:P,kk).*conj(sig_bb(1:P,kk)));
                elseif abs(kk-ll)<OverSampling
                    CFO = exp(1j*((kk-ll)/OverSampling*2*pi*(1:P)'));
                    % Assuming band spacing DeltaF, such that fs/DeltaF =
                    % OverSampling
                    sig_mix_u(d+1:d+P) = sig_mix_u(d+1:d+P)...
                        + sqrt(Pc*sig_pow(kk))*CAL_u(1:P,ll).*real(sig_bb(1:P,kk).*CFO * Phi(kk))...
                        + sqrt(sig_pow(kk)*sig_pow(ll))*real(sig_bb(1:P,kk).*conj(sig_bb(1:P,ll)).*CFO * Phi(kk) * conj(Phi(ll)));
                    sig_mix_l(d+1:d+P) = sig_mix_l(d+1:d+P)...
                        + sqrt(Pc*sig_pow(kk))*CAL_l(1:P,ll).*real(sig_bb(1:P,kk).*CFO * Phi(kk))...
                        + sqrt(sig_pow(kk)*sig_pow(ll))*real(sig_bb(1:P,kk).*conj(sig_bb(1:P,ll)).*CFO * Phi(kk) * conj(Phi(ll)));
                end
            end
        end      
        
        % Digital multiplication part, delay offset is adjusted
        bhi = fir1(200,1.1/OverSampling,'low');
        for bb=1:blk_num
            yout_u = filter(bhi,1,sqrt(2/Pc)*sig_mix_u(dOS+1:dOS+P).*CAL_u(1:P,bb));
            yout_l = filter(bhi,1,sqrt(2/Pc)*sig_mix_l(dOS+1:dOS+P).*CAL_l(1:P,bb));
            for nn = 1:length(Nrange)
                Sam_Num = Nrange(nn);
                results(onoffindex,bb,MCindex,nn) = yout_u(201:(200+Sam_Num))'*...
                    yout_l(201:(200+Sam_Num))/Sam_Num;
            end
        end
        
        % Double-arm multiplication for tagging
        for bb=1:blk_num
            for nn = 1:length(Nrange)
                Sam_Num = Nrange(nn);
                sumpow_hat = sum(squeeze(results(onoffindex,:,MCindex,nn)));
                sumpow_hat_res(onoffindex,bb,MCindex,nn) = sumpow_hat;
                sigma_hat = sumpow_hat*(blk_num-3)/OverSampling...
                                /sqrt(Sam_Num/OverSampling);
                TH = qfuncinv(0.05)*sigma_hat;
                ptag(onoffindex,bb,MCindex,nn) = results(onoffindex,bb,MCindex,nn)>TH;
            end
        end
    end
end

%% Plotting tagging rate curve
figure
for onoffindex = 1:blk_num+1
subplot(2,4,onoffindex)
for bb=1:blk_num
    plot(Nrange/(6*OverSampling),mean(squeeze(ptag(onoffindex,bb,:,:)),1),'o-','linewidth',3,'markersize',4);hold on
end
plot(Nrange/(6*OverSampling),ones(length(Nrange),1)*0.05,'k--','linewidth',3);hold on
ylim([0,1]);
grid on
title(['J' num2str(onoffindex) ' is off'])
xlabel('Time (\mus)')
ylabel('Probability')
end
legend('J1','J2','J3','J4','J5','J6')
%% test distribution of ED1
% onoffindex = 3;
% nn=6;
% Sam_Num = Nrange(nn);
% % Sam_Num = M;
% sig_pow_dB = [-5,-10,-15,-20,-25,-30]';
% sig_pow = 10.^(sig_pow_dB/10);
% sig_pow0 = sig_pow1;
% sig_pow0(onoffindex) = 0;
%     
% s_pow = sig_pow1(onoffindex);
% n_pow = sum(sig_pow0)/16*(6-1);
% N_prime = Sam_Num/SampleNumberAve(1);
% 
% % H0 emprical
% [f,xi] = ksdensity(squeeze(results_H0(onoffindex,:,nn)));
% mu = mean(squeeze(results_H0(onoffindex,:,nn)));
% sigma = sqrt(var(squeeze(results_H0(onoffindex,:,nn))));
% pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
% figure(99)
% % plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
% plot(xi,f,'--','linewidth',2);hold on
% 
% % H0 theoretical
% mu = 0;
% sigma = sqrt(n_pow^2/(1*N_prime));
% pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
% figure(99)
% plot(pdfxdata, normpdf(pdfxdata,mu,sigma),'linewidth',2);hold on
% 
% % H1 emprical
% [f,xi] = ksdensity(squeeze(results(onoffindex,:,nn)));
% mu = mean(squeeze(results(onoffindex,:,nn)));
% sigma = sqrt(var(squeeze(results(onoffindex,:,nn))));
% pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
% figure(99)
% % plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
% plot(xi,f,'--','linewidth',2)
% 
% 
% % H1 theoretical
% mu = s_pow;
% sigma = sqrt((2*s_pow^2+2*s_pow*n_pow+n_pow^2)/(1*N_prime));
% pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
% figure(99)
% plot(pdfxdata, normpdf(pdfxdata,mu,sigma),'linewidth',2);hold on
% grid on

%
% SINR_test = s_pow/n_pow;
% TH_test = qfuncinv(0.1)*n_pow/sqrt(N_prime);
% pd_test = qfunc((qfuncinv(0.1)-sqrt(N_prime)*SINR_test)/(sqrt((2*SINR_test^2+2*SINR_test+1))));
%

%% theoretical SINR verification
% Sam_Num = Sam_Num;
% SampleNumberAve_ED0 = 32;
% SampleNumberAve_ED1 = 32;
% SampleNumberAve_ED2 = 32;
% N_prime_ED0 = Sam_Num/(SampleNumberAve_ED0*1.1);
% N_prime_ED1 = Sam_Num/SampleNumberAve_ED1;
% N_prime_ED2 = Sam_Num/(SampleNumberAve_ED0*1.1);
% 
% blk_num = 6;
% pow_pool_fine = -10-linspace(0,20,100);
% % pow_pool_fine = -18;
% for pow_index = 1: 100
%     s_pow_ED1(pow_index) = 10^(pow_pool_fine(pow_index)/10)*1;
%     n_pow_ED1 = 1/SampleNumberAve_ED1*(blk_num-1)*1.1;
%     SINR_ED1(pow_index) = s_pow_ED1(pow_index)/n_pow_ED1;
%     pd_theo_ED1_new(pow_index) = qfunc((qfuncinv(0.1)-sqrt(N_prime_ED1)*SINR_ED1(pow_index))/(sqrt(1+2*SINR_ED1(pow_index)+2*SINR_ED1(pow_index)^2)));
% end
% 
% figure(1)
% plot(pow_pool_fine,pd_theo_ED1_new,'linewidth',3);hold on

