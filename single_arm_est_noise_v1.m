%  08/15/2016
%  estimate interference floor and use it to set threshold for identifying
%  dominant blockers

clear;clc;%clf;close all


% rand('seed',3)
blk_num = 6;
M = 1e2;
N = 63;
SampleNumberAve = 32;
stat_num = 2e2;
P = M*N;
noise_est_ratio = 2;
P_noise_est = P*noise_est_ratio;
L = M*SampleNumberAve*stat_num;
Sam_Num = 5e3;
Sam_Num_noise_est = noise_est_ratio*Sam_Num;

for ii = 1:blk_num
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
%%
other_pow = 1e-4;
sig_pow0 = ones(blk_num,1)*other_pow;
sig_pow1 = ones(blk_num,1)*other_pow;

pow_pool = -10-linspace(0,18,10);
% pow_pool = -15;
pow_num = 10;
TH = zeros(pow_num,1);
pd = zeros(pow_num,1);

for pow_index = 1:pow_num
j1_dB = pow_pool(pow_index);
sig_pow1(4) = 1;
sig_pow0(4) = 1;
% sig_pow1(2:6) = ones(5,1)*0.1;
% sig_pow0(2:6) = ones(5,1)*0.1;
sig_pow1(1) = 10^(j1_dB/10);
sig_pow0(1) = 0;

%% calibration sequences
% cal0=PNgenerator_v5(N,N,1);
% CAL_noperm=LowerRate_v2(cal0,P);
% CAL=zeros(P,N);
% for mm=1:M
%     CAL((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
% end
%% noise power estimation
runtimesTH = 500;
b_index_est = 2;
for runindex = 1:runtimesTH
    CAL = randi(2,P_noise_est,blk_num)*2-3;
    CAL_DSP = CAL(:,b_index_est);
    CAL(:,b_index_est) = 0;
    sig_cache = zeros(P_noise_est,blk_num);
    start_point = randi(L-P_noise_est-1,blk_num,1);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii,1):start_point(ii,1)+P_noise_est-1,randi(blk_num));
    end
    sig_mix0_part1 = zeros(P_noise_est,1);
    sig_mix1_part1 = zeros(P_noise_est,1);

    sig_mix0_part2 = zeros(P_noise_est,1);
    sig_mix1_part2 = zeros(P_noise_est,1);

    sig_mix0 = zeros(P_noise_est,1);
    sig_mix1 = zeros(P_noise_est,1);
    PhaseShift = exp(1j*rand(1,blk_num));
    
    for ii=1:blk_num
        sig_mix0_part1 = sig_mix0_part1+sig_cache(:,ii).*CAL(:,ii)*PhaseShift(ii)*sqrt(sig_pow0(ii));
        sig_mix1_part1 = sig_mix1_part1+sig_cache(:,ii).*CAL(:,ii)*PhaseShift(ii)*sqrt(sig_pow1(ii));
    end
    
    for ii=1:blk_num
        for jj=1:blk_num
            if ii ~= jj
            carrier = exp(1j*((ii-jj)/32*2*pi*(1:P_noise_est)'))*PhaseShift(ii);
            sig_mix0_part2 = sig_mix0_part2 + sqrt(sig_pow0(ii))*CAL(:,jj).*...
                real(sig_cache(:,ii).*carrier);
            sig_mix1_part2 = sig_mix1_part2 + sqrt(sig_pow1(ii))*CAL(:,jj).*...
                real(sig_cache(:,ii).*carrier);
            end
        end
    end
    sig_mix0 = sig_mix0_part1+sig_mix0_part2;
    sig_mix1 = sig_mix1_part1+sig_mix1_part2;
    
    
    bhi = fir1(100,1/SampleNumberAve,'low');
    sigma_hat1(pow_index,runindex) = norm(filter(bhi,1,sig_mix1(1:Sam_Num_noise_est)...
            .*CAL_DSP(1:Sam_Num_noise_est)))/sqrt(Sam_Num_noise_est);
    sigma_hat0(pow_index,runindex) = norm(filter(bhi,1,sig_mix0(1:Sam_Num_noise_est)...
            .*CAL_DSP(1:Sam_Num_noise_est)))/sqrt(Sam_Num_noise_est);
end

%% tagging
runtimes = 5e2;
for runindex=1:runtimes
%     CAL=zeros(P,N);
%     for mm=1:M
%         CAL((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
%     end
    CAL = randi(2,P,blk_num)*2-3;

    %CAL=CAL(:,randperm(N));
    sig_cache = zeros(P,blk_num);
    start_point = randi(L-P-1,blk_num,1);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii,1):start_point(ii,1)+P-1,randi(blk_num));
    end

    sig_mix0_part1 = zeros(P,1);
    sig_mix1_part1 = zeros(P,1);

    sig_mix0_part2 = zeros(P,1);
    sig_mix1_part2 = zeros(P,1);

    sig_mix0 = zeros(P,1);
    sig_mix1 = zeros(P,1);
    PhaseShift = exp(1j*rand(1,blk_num));
    
    
    for ii=1:blk_num
        sig_mix0_part1 = sig_mix0_part1+sig_cache(:,ii).*CAL(:,ii)*PhaseShift(ii)*sqrt(sig_pow0(ii));
        sig_mix1_part1 = sig_mix1_part1+sig_cache(:,ii).*CAL(:,ii)*PhaseShift(ii)*sqrt(sig_pow1(ii));
    end

    for ii=1:blk_num
        for jj=1:blk_num
            if ii ~= jj
            carrier = exp(1j*((ii-jj)/32*2*pi*(1:P)'))*PhaseShift(ii);
            sig_mix0_part2 = sig_mix0_part2 + sqrt(sig_pow0(ii))*CAL(:,jj).*...
                real(sig_cache(:,ii).*carrier);
            sig_mix1_part2 = sig_mix1_part2 + sqrt(sig_pow1(ii))*CAL(:,jj).*...
                real(sig_cache(:,ii).*carrier);
            end
        end
    end
    sig_mix0 = sig_mix0_part1+sig_mix0_part2;
    sig_mix1 = sig_mix1_part1+sig_mix1_part2;
    
    %% original algorithm
        
        bhi = fir1(100,1/SampleNumberAve,'low');

        ED_results1(pow_index,runindex) = norm(filter(bhi,1,sig_mix1(1:Sam_Num)...
            .*CAL(1:Sam_Num,1)))/sqrt(Sam_Num);
        ED_results0(pow_index,runindex) = norm(filter(bhi,1,sig_mix0(1:Sam_Num)...
            .*CAL(1:Sam_Num,1)))/sqrt(Sam_Num);
     
    
end
%
temp = sort(ED_results0(pow_index,:),'ascend');
TH_opt(pow_index) = temp(runtimes*0.9);
pd_opt(pow_index) = sum(ED_results1(pow_index,:)>TH_opt(pow_index))/runtimes;


% temp = sort(ED_results2(1,:),'ascend');
% TH2(pow_index) = temp(runtimes*0.9);
% pd2(pow_index) = sum(ED_results2(2,:)>TH2(pow_index))/runtimes;
end
%%
for pow_index = 1:pow_num
    for runindexTH = 1:runtimesTH
%         pd_noise_est_temp(runindexTH) = sum(ED_results1(pow_index,:)...
%             >(1+qfuncinv(0.1)*sqrt(2*8/Sam_Num))*sigma_hat1(pow_index,runindexTH))/runtimes;
%         pfa_noise_est_temp(runindexTH) = sum(ED_results0(pow_index,:)...
%             >(1+qfuncinv(0.1)*sqrt(2*8/Sam_Num))*sigma_hat0(pow_index,runindexTH))/runtimes;
        pd_noise_est_temp(runindexTH) = sum(ED_results1(pow_index,:)...
            >1.085*sigma_hat1(pow_index,runindexTH))/runtimes;
        pfa_noise_est_temp(runindexTH) = sum(ED_results0(pow_index,:)...
            >1.085*sigma_hat0(pow_index,runindexTH))/runtimes;
    
    end
    pd_power_est(pow_index) = sum(pd_noise_est_temp)/runtimesTH;
    pfa_power_est(pow_index) = sum(pfa_noise_est_temp)/runtimesTH;
end

%
figure
plot(pow_pool,pd_opt);hold on
plot(pow_pool,pd_power_est);hold on
plot(pow_pool,pfa_power_est);hold on
grid on
legend('PD, optimal TH','PD, estimated TH','PFA,estimated TH')
ylim([0,1])
%%
ii=5;
[f,xi] = ksdensity(ED_results1(ii,:));
mu = mean(ED_results1(ii,:))
sigma = sqrt(var(ED_results1(ii,:)))
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure
% plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
plot(xi,f);hold on

[f,xi] = ksdensity(ED_results0(ii,:));
mu = mean(ED_results0(ii,:))
sigma = sqrt(var(ED_results0(ii,:)))
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
% plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
plot(xi,f)

mu = mean(sigma_hat0(ii,:))
sigma = sqrt(var(sigma_hat0(ii,:)))
[f,xi] = ksdensity(sigma_hat0(ii,:));
plot(xi,f,'--')

mu = mean(sigma_hat1(ii,:))
sigma = sqrt(var(sigma_hat1(ii,:)))
[f,xi] = ksdensity(sigma_hat1(ii,:));
plot(xi,f,'--')
legend('Meaasured Power (H1)','Meaasured Power (H0)','Est. Noise Pow (H0)','Est. Noise Pow (H1)')

%%
% figure
% subplot(211)
% plot(sort(ED_results1(1,:)));hold on
% plot(sort(ED_results1(2,:)));
% subplot(212)
% plot(sort(ED_results2(1,:)));hold on
% plot(sort(ED_results2(2,:)));
%%

% figure
% subplot(211)
% plot(abs(temp(51:1050)));hold on
% plot(abs(sig_cache(1:1000,4)))
% subplot(212)
% plot(abs(temp(51:1050)-sig_cache(1:1000,4)*PhaseShift(4)));
%%
% figure
% subplot(211)
% plot(angle(temp(51:1050)));hold on
% plot(angle(sig_cache(1:1000,4)*PhaseShift(4)))
%%
% figure
% subplot(311)
% plot(abs(sig_mix0(1:5000)));hold on
% subplot(312)
% plot(abs(y_j4_hat));
% subplot(313)
% plot(abs(y_j4_hat-sig_mix0(1:5000)));
