%  01/15/2017
%  There are 6 cases. In i^th case blocker i is turned on and off to test
%  Pd and Pfa rate. The proposed TH is tested which provides less than 0.05
%  flase alarm rate but less pd as well. Optimal TH is used to generate
%  0.05 pfa and its pd serves as benchmark

clear;clc;%clf;close all
rng('default')
blk_num = 6;
% rand('seed',3)
% oversampling ratio, which is PN rate v.s. blocker rate
OverSampling = 8;

total_time = 1e-3; % maximum running time is 1ms
Lmax = fix(total_time/(1/6e6/OverSampling));

% Length of PN per period
M = 1600;
N = fix(Lmax/M);
P = M*N;
L = 100*P;

Nrange = fix(P*(0.1:0.1:0.9));
Pc = 10;

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

% pow_num = 11;
% pow_pool = -14;
% pow_pool = -10-linspace(0,20,pow_num);
% TH = zeros(pow_num,1);
% pd = zeros(pow_num,1);
%%
runtimes = 5e2;
results_H0 = zeros(blk_num,runtimes,length(Nrange));
results_H1 = zeros(blk_num,runtimes,length(Nrange));
pfa = zeros(1,length(Nrange));
pd = zeros(1,length(Nrange));
pd_opt = zeros(1,length(Nrange));

for onoffindex = 1:blk_num

    XXX = ['simulating, ',num2str(onoffindex/blk_num,2),' finished'];
    disp(XXX)
    
    % blocker power settings
    sig_pow_dB = [-3,-6,-9,-12,-15,-18]';
    
    sig_pow1 = 10.^(sig_pow_dB/10);
    sig_pow0 = sig_pow1;
    sig_pow0(onoffindex) = 0;
    
    % tagging
    
    for runindex=1:runtimes

%         CAL_u = randi(2,P,blk_num)*2-3;
%         CAL_l = randi(2,P,blk_num)*2-3;
        CALseed = randi(2,M+1,blk_num)*2-3;
        CAL_u0 = kron(ones(N,1),CALseed);
        CAL_u = CAL_u0(1:P,:);
        CALseed = randi(2,M+2,blk_num)*2-3;
        CAL_l0 = kron(ones(N,1),CALseed);
        CAL_l = CAL_l0(1:P,:);

        start_point = randi(L-P-1,blk_num,1);
        for ii=1:blk_num
            sig_bb(:,ii) = sig(start_point(ii,1):start_point(ii,1)+P-1,randi(blk_num));
        end

        sig_mixH0_u = zeros(2*P,1);
        sig_mixH1_u = zeros(2*P,1);
        sig_mixH0_l = zeros(2*P,1);
        sig_mixH1_l = zeros(2*P,1);

        Phi = exp(1j*rand(1,blk_num));
        
        % pilots contains two CAL passing through diode
%         d = randi(5e3);
%         pilot_u = 0.5*Pc*[zeros(d,1);CALseed(:,1)].*conj([zeros(d,1);CALseed(:,1)])...
%             +0.5*Pc*[zeros(d,1);CALseed(:,2)].*conj([zeros(d,1);CALseed(:,2)])...
%             +Pc*real([zeros(d,1);CALseed(:,1)].*conj([zeros(d,1);CALseed(:,2)]));
%         [peak,dhat] = max(abs(conv(pilot_u,flipud(CALseed(:,1).*CALseed(:,2)))));
%         dOS = dhat-size(CALseed(:,1),1);
        d = 0;
        dOS = 0;
        

        for kk=1:blk_num
            for ll=1:blk_num
                if kk==ll
                    % Hypothesis 0 part
                    sig_mixH0_u(d+1:d+P) = sig_mixH0_u(d+1:d+P)...
                        +sqrt(Pc*sig_pow0(kk))*real(sig_bb(1:P,kk)*Phi(kk)).*CAL_u(1:P,kk)...
                        +0.5*sig_pow0(kk)*(sig_bb(1:P,kk).*conj(sig_bb(1:P,kk)));
                    sig_mixH0_l(d+1:d+P) = sig_mixH0_l(d+1:d+P)...
                        +sqrt(Pc*sig_pow0(kk))*real(sig_bb(1:P,kk)*Phi(kk)).*CAL_l(1:P,kk)...
                        +0.5*sig_pow0(kk)*(sig_bb(1:P,kk).*conj(sig_bb(1:P,kk)));
                    % Hypothesis 1 part
                    sig_mixH1_u(d+1:d+P) = sig_mixH1_u(d+1:d+P)...
                        +sqrt(Pc*sig_pow1(kk))*real(sig_bb(1:P,kk)*Phi(kk)).*CAL_u(1:P,kk)...
                        +0.5*sig_pow1(kk)*(sig_bb(1:P,kk).*conj(sig_bb(1:P,kk)));
                    sig_mixH1_l(d+1:d+P) = sig_mixH1_l(d+1:d+P)...
                        +sqrt(Pc*sig_pow1(kk))*real(sig_bb(1:P,kk)*Phi(kk)).*CAL_l(1:P,kk)...
                        +0.5*sig_pow1(kk)*(sig_bb(1:P,kk).*conj(sig_bb(1:P,kk)));
                elseif abs(kk-ll)<OverSampling
                    CFO = exp(1j*((kk-ll)/OverSampling*2*pi*(1:P)'));
                    % Assuming band spacing DeltaF, such that fs/DeltaF = 32
                    % H0 part
                    sig_mixH0_u(d+1:d+P) = sig_mixH0_u(d+1:d+P)...
                        + sqrt(Pc*sig_pow0(kk))*CAL_u(1:P,ll).*real(sig_bb(1:P,kk).*CFO * Phi(kk))...
                        + sqrt(sig_pow0(kk)*sig_pow0(ll))*real(sig_bb(1:P,kk).*conj(sig_bb(1:P,ll)).*CFO * Phi(kk) * conj(Phi(ll)));
                    sig_mixH0_l(d+1:d+P) = sig_mixH0_l(d+1:d+P)...
                        + sqrt(Pc*sig_pow0(kk))*CAL_l(1:P,ll).*real(sig_bb(1:P,kk).*CFO * Phi(kk))...
                        + sqrt(sig_pow0(kk)*sig_pow0(ll))*real(sig_bb(1:P,kk).*conj(sig_bb(1:P,ll)).*CFO * Phi(kk) * conj(Phi(ll)));
                    % H1 part
                    sig_mixH1_u(d+1:d+P) = sig_mixH1_u(d+1:d+P)...
                        + sqrt(Pc*sig_pow1(kk))*CAL_u(1:P,ll).*real(sig_bb(1:P,kk).*CFO * Phi(kk))...
                        + sqrt(sig_pow1(kk)*sig_pow1(ll))*real(sig_bb(1:P,kk).*conj(sig_bb(1:P,ll)).*CFO * Phi(kk) * conj(Phi(ll)));
                    sig_mixH1_l(d+1:d+P) = sig_mixH1_l(d+1:d+P)...
                        + sqrt(Pc*sig_pow1(kk))*CAL_l(1:P,ll).*real(sig_bb(1:P,kk).*CFO * Phi(kk))...
                        + sqrt(sig_pow1(kk)*sig_pow1(ll))*real(sig_bb(1:P,kk).*conj(sig_bb(1:P,ll)).*CFO * Phi(kk) * conj(Phi(ll)));
                end
            end
        end      
        
        bhi = fir1(200,1.1/OverSampling,'low');
        yout_u1 = filter(bhi,1,sqrt(2/Pc)*sig_mixH1_u(dOS+1:dOS+P).*CAL_u(1:P,onoffindex));
        yout_l1 = filter(bhi,1,sqrt(2/Pc)*sig_mixH1_l(dOS+1:dOS+P).*CAL_l(1:P,onoffindex));
        for nn = 1:length(Nrange)
            Sam_Num = Nrange(nn);
            results_H1(onoffindex,runindex,nn) = yout_u1(201:(200+Sam_Num))'*...
                yout_l1(201:(200+Sam_Num))/Sam_Num;
        end

        yout_u0 = filter(bhi,1,sqrt(2/Pc)*sig_mixH0_u(dOS+1:dOS+P).*CAL_u(1:P,onoffindex));
        yout_l0 = filter(bhi,1,sqrt(2/Pc)*sig_mixH0_l(dOS+1:dOS+P).*CAL_l(1:P,onoffindex));
        for nn = 1:length(Nrange)
            Sam_Num = Nrange(nn);
            results_H0(onoffindex,runindex,nn) = yout_u0(201:(200+Sam_Num))'*...
                yout_l0(201:(200+Sam_Num))/Sam_Num;
        end
    end


    for nn = 1:length(Nrange)
        Sam_Num = Nrange(nn);
%         sigma_hat(onoffindex,nn) = 1*(blk_num-1)/OverSampling...
%                         /sqrt(Sam_Num/OverSampling);
        sigma_hat(nn) = sum(sig_pow1)*(blk_num-1)/OverSampling...
                        /sqrt(Sam_Num/OverSampling);
        TH(nn) = qfuncinv(0.05)*sigma_hat(nn);
        pfa(onoffindex,nn) = sum(squeeze(results_H0(onoffindex,:,nn))...
                                    >TH(nn))/runtimes;
        pd(onoffindex,nn) = sum(squeeze(results_H1(onoffindex,:,nn))...
                                    >TH(nn))/runtimes;

        temp = sort(squeeze(results_H0(onoffindex,:,nn)),'ascend');
        TH_opt = temp(runtimes*0.95);
        pd_opt(onoffindex,nn) = sum(squeeze(results_H1(onoffindex,:,nn))>TH_opt)/runtimes;

    end

end

%% Plotting Pd & Pfa curve
figure
for onoffindex = 1:blk_num
subplot(2,3,onoffindex)
plot(Nrange/(6*OverSampling),pd(onoffindex,:),'o-','linewidth',3,'markersize',4);hold on
plot(Nrange/(6*OverSampling),pfa(onoffindex,:),'o-','linewidth',3,'markersize',4);hold on
plot(Nrange/(6*OverSampling),pd_opt(onoffindex,:),'o--','linewidth',3,'markersize',4);hold on
plot(Nrange/(6*OverSampling),ones(length(Nrange),1)*0.05,'k--','linewidth',3);hold on
ylim([0,1]);
grid on
title(['J' num2str(onoffindex) ' with ' num2str(sig_pow_dB(onoffindex)) 'dBm'])
xlabel('Time (\mus)')
ylabel('Probability')
end
legend('PD (proposed TH)','PFA (proposed TH)','PD (opt TH)','PFA Target')

%% test distribution of ED1
onoffindex = 5;
nn=3;
Sam_Num = Nrange(nn);
% Sam_Num = M;
% sig_pow_dB = [-3,-6,-15,-20,-25,-30]';
sig_pow1 = 10.^(sig_pow_dB/10);
sig_pow0 = sig_pow1;
sig_pow0(onoffindex) = 0;
    
s_pow = sig_pow1(onoffindex);
n_pow = sum(sig_pow0)/OverSampling*(6-1);
% N_prime = Sam_Num/OverSampling;
N_prime = Sam_Num/OverSampling/2;

% H0 emprical
[f,xi] = ksdensity(squeeze(results_H0(onoffindex,:,nn)));
mu = mean(squeeze(results_H0(onoffindex,:,nn)));
sigma = sqrt(var(squeeze(results_H0(onoffindex,:,nn))));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(99)
% plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
plot(xi,f,'--','linewidth',2);hold on

% H0 theoretical
mu = 0;
sigma = sqrt(n_pow^2/(1*N_prime));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(99)
plot(pdfxdata, normpdf(pdfxdata,mu,sigma),'linewidth',2);hold on

% H1 emprical
[f,xi] = ksdensity(squeeze(results_H1(onoffindex,:,nn)));
mu = mean(squeeze(results_H1(onoffindex,:,nn)));
sigma = sqrt(var(squeeze(results_H1(onoffindex,:,nn))));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(99)
% plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
plot(xi,f,'--','linewidth',2)


% H1 theoretical
mu = s_pow;
sigma = sqrt((2*s_pow^2+2*s_pow*n_pow+n_pow^2)/(1*N_prime));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(99)
plot(pdfxdata, normpdf(pdfxdata,mu,sigma),'linewidth',2);hold on
grid on

%
% SINR_test = s_pow/n_pow;
% TH_test = qfuncinv(0.1)*n_pow/sqrt(N_prime);
% pd_test = qfunc((qfuncinv(0.1)-sqrt(N_prime)*SINR_test)/(sqrt((2*SINR_test^2+2*SINR_test+1))));
%

%% theoretical SINR verification
% Sam_Num = Sam_Num;
% OverSampling_ED0 = 16;
% OverSampling_ED1 = 16;
% OverSampling_ED2 = 16;
% N_prime_ED0 = Sam_Num/(OverSampling_ED0*1.1);
% N_prime_ED1 = Sam_Num/OverSampling_ED1;
% N_prime_ED2 = Sam_Num/(OverSampling_ED0*1.1);
figure
blk_num = 6;
Nrange_fine = linspace(Nrange(1),Nrange(end),100);

for onoffindex = 1:blk_num
    s_pow = 10^(sig_pow_dB(onoffindex)/10);
    sig_pow1 = 10.^(sig_pow_dB/10);
    sig_pow0 = sig_pow1;
    sig_pow0(onoffindex) = 0;
    n_pow = sum(sig_pow0)/OverSampling*(blk_num-1);
    for n_index = 1: 100
        N_prime = Nrange_fine(n_index)/OverSampling/2;
        SINR(n_index) = s_pow/n_pow;
        pd_theo(n_index) = qfunc((qfuncinv(0.05)-sqrt(N_prime)*SINR(n_index))/(sqrt(1+2*SINR(n_index)+2*SINR(n_index)^2)));
    end
    subplot(2,3,onoffindex)
    plot(Nrange_fine/(6*OverSampling),pd_theo,'linewidth',3);hold on
    ylim([0,1])
end
