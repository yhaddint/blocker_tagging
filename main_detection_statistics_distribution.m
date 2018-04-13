%  08/15/2016
%  estimate interference floor and use it to set threshold for identifying
%  dominant blockers
% It verifies theoretical analysis on distribution of detection statistics 

clear;clc;%clf;close all


% rand('seed',3)
blk_num = 6;
M = 1e2;
N = 63;
SampleNumberAve = 32;
stat_num = 2e2;
P = M*N;
L = M*SampleNumberAve*stat_num;
pow_PN = 10^(50/10);

for kk = 1:blk_num
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
    sig(:,kk) = temp_data(end-L+1:end)./sqrt(temp_data(end-L+1:end)'*temp_data(end-L+1:end)/L);
end

pow_pool = -10-linspace(0,18,10);
% pow_pool = -15
pow_num = 10;
TH = zeros(pow_num,1);
pd = zeros(pow_num,1);

for pow_index = 1:pow_num
j1_dB = pow_pool(pow_index);
sig_pow0 = zeros(blk_num,1);
sig_pow1 = zeros(blk_num,1);

% power of other interferers
% sig_pow1(4) = 1;
% sig_pow0(4) = 1;
sig_pow1(2:6) = ones(5,1)*0.2;
sig_pow0(2:6) = ones(5,1)*0.2;
sig_pow1(1) = 10^(j1_dB/10);
sig_pow0(1) = 0;

%% calibration sequences
% cal0=PNgenerator_v5(N,N,1);
% CAL_noperm=LowerRate_v2(cal0,P);
% CAL=zeros(P,N);
% for mm=1:M
%     CAL((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
% end

%% tagging
runtimes = 5e2;
for runindex=1:runtimes
%     CAL=zeros(P,N);
%     for mm=1:M
%         CAL((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
%     end
    CAL = randi(2,P,blk_num)*2-3;

    %CAL=CAL(:,randperm(N));
    start_point = randi(L-P-1,blk_num,1);
    for kk=1:blk_num
        sig_cache(:,kk) = sig(start_point(kk,1):start_point(kk,1)+P-1,randi(blk_num));
    end
%     sig_cache(:,4) = sig_other(start_point(ii,1):start_point(ii,1)+P-1,1);
    
    sig_mix0 = zeros(P,1);
    sig_mix1 = zeros(P,1);

    PhaseShift = exp(1j*rand(1,blk_num));

    for kk=1:blk_num
        for ll=1:kk
            if kk==ll
                sig_mix0 = sig_mix0...
                    +sqrt(sig_pow0(kk))*real(sig_cache(:,kk)*PhaseShift(kk)).*CAL(:,kk)...
                    +0.5*sig_pow0(kk)/sqrt(pow_PN)*conj(sig_cache(:,kk)).*sig_cache(:,kk);
                sig_mix1 = sig_mix1...
                    +sqrt(sig_pow1(kk))*real(sig_cache(:,kk)*PhaseShift(kk)).*CAL(:,kk)...
                    +0.5*sig_pow1(kk)/sqrt(pow_PN)*conj(sig_cache(:,kk)).*sig_cache(:,kk);
            else
                carrier = exp(1j*((kk-ll)/SampleNumberAve*2*pi*(1:P)'))*PhaseShift(kk);
                sig_mix0 = sig_mix0...
                    +sqrt(sig_pow0(kk))*CAL(:,ll).*real(sig_cache(:,kk).*carrier)...
                    +sqrt(sig_pow0(kk)*sig_pow0(ll)/pow_PN)*real(sig_cache(:,kk).*conj(sig_cache(:,ll)).*carrier);
                sig_mix1 = sig_mix1...
                    +sqrt(sig_pow1(kk))*CAL(:,ll).*real(sig_cache(:,kk).*carrier)...
                    +sqrt(sig_pow1(kk)*sig_pow1(ll)/pow_PN)*real(sig_cache(:,kk).*conj(sig_cache(:,ll)).*carrier);
            end
        end
    end
    
    %% original algorithm
    Sam_Num = 5e3;
    bhi = fir1(200,1/SampleNumberAve+1/64,'low');
    LPF_out1 = filter(bhi,1,sig_mix1.*CAL(:,1));
    LPF_out0 = filter(bhi,1,sig_mix0.*CAL(:,1));
    ED_results1(pow_index,runindex) = norm(LPF_out1(201:(200+Sam_Num)))^2/Sam_Num;
    ED_results0(pow_index,runindex) = norm(LPF_out0(201:(200+Sam_Num)))^2/Sam_Num; 

end
%
temp = sort(ED_results0(pow_index,:),'ascend');
TH(pow_index) = temp(runtimes*0.9);
pd(pow_index) = sum(ED_results1(pow_index,:)>TH(pow_index))/runtimes;

% temp = sort(ED_results2(1,:),'ascend');
% TH2(pow_index) = temp(runtimes*0.9);
% pd2(pow_index) = sum(ED_results2(2,:)>TH2(pow_index))/runtimes;
% if pow_index==1
%    1 
% end
end
%%

figure(1)
plot(pow_pool,pd,'--o','linewidth',3);hold on
% plot(pow_pool,pd2);hold on
grid on
% legend('algorithm 1','algorithm 2')
%%
% pow_index = 1;
% figure
% plot(sort(ED_results0(pow_index,:)));hold on
% plot(sort(ED_results1(pow_index,:)));
%%
pow_index = 6;
[f,xi] = ksdensity(ED_results0(pow_index,:));
mu = mean(ED_results0(pow_index,:));
sigma = sqrt(var(ED_results0(pow_index,:)));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure
plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
plot(xi,f)

[f,xi] = ksdensity(ED_results1(pow_index,:));
mu = mean(ED_results1(pow_index,:));
sigma = sqrt(var(ED_results1(pow_index,:)));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
plot(xi,f)
%%
% ydata = abs(fft(sig_mix1(1:Sam_Num).*CAL(1:Sam_Num,1)));
% ydata = abs(fft(sig_cache(:,1)));
% xdata = linspace(-pi,pi,Sam_Num);
% figure
% plot(xdata,[ydata(2501:5000);ydata(1:2500)])
%% theoretical SINR verification
Sam_Num = 5e3;
SampleNumberAve = 32;
N_prime = Sam_Num/SampleNumberAve;
W_ratio = SampleNumberAve;
blk_num = 6;
pow_pool_fine = -10-linspace(0,18,100);
for pow_index = 1: 100
    s_pow(pow_index) = 10^(pow_pool_fine(pow_index)/10)*(1/2+1/2/W_ratio*(blk_num-1));
    n_pow = 1/2/W_ratio*(blk_num-1);
    SINR(pow_index) = s_pow(pow_index)/n_pow;
    pd_theo_ED0(pow_index) = qfunc((qfuncinv(0.1)-sqrt(Sam_Num/2)*SINR(pow_index))/(1+SINR(pow_index)));
    
    TH_theo_ED1(pow_index) = qfuncinv(0.1)*n_pow/sqrt(N_prime);
    pd_theo_ED1(pow_index) = qfunc((-s_pow(pow_index)+TH_theo_ED1(pow_index))/(sqrt((2*s_pow(pow_index)^2+2*s_pow(pow_index)*n_pow+n_pow^2)/N_prime)));
    
    pd_theo_ED2(pow_index) = qfunc((qfuncinv(0.1)-sqrt(N_prime/2)*SINR(pow_index))/(1+2*SINR(pow_index)));
    pd_theo_ED0_new(pow_index) = qfunc((qfuncinv(0.1)-sqrt(N_prime/2)*SINR(pow_index))/(1+SINR(pow_index)));
    pd_theo_ED1_new(pow_index) = qfunc((qfuncinv(0.1)-sqrt(N_prime)*SINR(pow_index))/(sqrt(1+2*SINR(pow_index)+2*SINR(pow_index)^2)));
    pd_theo_ED2_new(pow_index) = qfunc((qfuncinv(0.1)-sqrt(2*N_prime)*SINR(pow_index))/(1+2*SINR(pow_index)));

end
% figure(1)
% plot(pow_pool_fine,pd_theo_ED0);hold on
% plot(pow_pool_fine,pd_theo_ED1);hold on
% plot(pow_pool_fine,pd_theo_ED2);hold on
% legend('ED0','ED1','ED2')
figure(1)
plot(pow_pool_fine,pd_theo_ED0_new,'linewidth',3);hold on
plot(pow_pool_fine,pd_theo_ED1_new,'linewidth',3);hold on
plot(pow_pool_fine,pd_theo_ED2_new,'linewidth',3);hold on
