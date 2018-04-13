%  08/15/2016
%  Performance analysis of three energy detection, 1) single-arm ED, 2)
%  double arm cross-multiplication ED, 3) optimal ED using double-arm
% Detection statistics are verified and Pd is evaluate from empical TH

clear;clc;%clf;close all


% rand('seed',3)
blk_num = 6;
M = 1e2;
N = 63;
SampleNumberAve = 32;
stat_num = 2e2;
P = M*N;
L = M*SampleNumberAve*stat_num;

% f = [0 1/32-1/64 1/32 1];
% mhi = [1 1 0 0];
% bhi = fir2(400,f,mhi);
bhi = fir1(400,1/SampleNumberAve*0.7,'low');
% freqz(bhi,1)

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
    temp1_data = filter(bhi,1,temp_data);
%     temp_data = data;
%     unnormalized = real(temp_data(end-L+1:end));
%     sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     unnormalized = imag(temp_data(end-L+1:end));
%     sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     sig_i(:,ii) = zeros(size(sig_r(:,ii)));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
    sig(:,ii) = temp1_data(end-L+1:end)./sqrt(temp1_data(end-L+1:end)'*temp1_data(end-L+1:end)/L);

end

pow_num = 11;
% pow_pool = -14;
pow_pool = -10-linspace(0,20,pow_num);
TH = zeros(pow_num,1);
pd = zeros(pow_num,1);
%%
for pow_index = 1:pow_num

XXX = ['simulating, ',num2str(pow_index/pow_num,2),' finished'];
disp(XXX)
    
j1_dB = pow_pool(pow_index);

sig_pow0 = zeros(blk_num,1);
sig_pow1 = zeros(blk_num,1);
sig_pow1(1) = 10^(j1_dB/10);

sig_pow1(4) = 1;
sig_pow0(4) = 1;

% sig_pow1(2:6) = ones(5,1)*0.2;
% sig_pow0(2:6) = ones(5,1)*0.2;

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
    CAL_u = randi(2,P,blk_num)*2-3;
%     CAL_u(:,[2,3,5,6]) = 0;
    CAL_l = randi(2,P,blk_num)*2-3;

    %CAL=CAL(:,randperm(N));
    start_point = randi(L-P-1,blk_num,1);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii,1):start_point(ii,1)+P-1,randi(blk_num));
    end

    sig_mix0_u = zeros(P,1);
    sig_mix1_u = zeros(P,1);
    sig_mix0_l = zeros(P,1);
    sig_mix1_l = zeros(P,1);
    
    PhaseShift = exp(1j*rand(1,blk_num));

    for kk=1:blk_num
        for ll=1:blk_num
            if kk==ll
                sig_mix0_u = sig_mix0_u+sqrt(2*sig_pow0(kk))*real(sig_cache(:,kk)*PhaseShift(kk)).*CAL_u(:,kk);
                sig_mix1_u = sig_mix1_u+sqrt(2*sig_pow1(kk))*real(sig_cache(:,kk)*PhaseShift(kk)).*CAL_u(:,kk);
                sig_mix0_l = sig_mix0_l+sqrt(2*sig_pow0(kk))*real(sig_cache(:,kk)*PhaseShift(kk)).*CAL_l(:,kk);
                sig_mix1_l = sig_mix1_l+sqrt(2*sig_pow1(kk))*real(sig_cache(:,kk)*PhaseShift(kk)).*CAL_l(:,kk);
            else
                carrier = exp(1j*((kk-ll)/SampleNumberAve*2*pi*(1:P)'))*PhaseShift(kk);
                % H0 part
                sig_mix0_u = sig_mix0_u + sqrt(2*sig_pow0(kk))*CAL_u(:,ll).*...
                    real(sig_cache(:,kk).*carrier);
                sig_mix0_l = sig_mix0_l + sqrt(2*sig_pow0(kk))*CAL_l(:,ll).*...
                    real(sig_cache(:,kk).*carrier);
                % H1 part
                sig_mix1_u = sig_mix1_u + sqrt(2*sig_pow1(kk))*CAL_u(:,ll).*...
                    real(sig_cache(:,kk).*carrier);
                sig_mix1_l = sig_mix1_l + sqrt(2*sig_pow1(kk))*CAL_l(:,ll).*...
                    real(sig_cache(:,kk).*carrier);
            end
        end
    end
    
    %% 3 different detector
    Sam_Num = 5e3;
    bhi = fir1(200,1.1/SampleNumberAve,'low');
%     freqz(bhi,1)
    yout_u = filter(bhi,1,sig_mix1_u.*CAL_u(:,1));
%     yout_u = sig_mix1_u.*CAL_u(:,1);

    yout_l = filter(bhi,1,sig_mix1_l.*CAL_l(:,1));
%     yout_l = sig_mix1_l.*CAL_l(:,1);
    results_SA_H1(pow_index,runindex) = yout_u(201:(200+Sam_Num))'*yout_u(201:(200+Sam_Num))/Sam_Num;
    results_DA_PED_H1(pow_index,runindex) = yout_u(201:(200+Sam_Num))'*yout_l(201:(200+Sam_Num))/Sam_Num;
    results_DA_OED_H1(pow_index,runindex) = norm((yout_u(201:(200+Sam_Num))+yout_l(201:(200+Sam_Num)))/2)^2/Sam_Num;

    yout_u = filter(bhi,1,sig_mix0_u.*CAL_u(:,1));
%     yout_u = sig_mix0_u.*CAL_u(:,1);
    yout_l = filter(bhi,1,sig_mix0_l.*CAL_l(:,1));
%     yout_l = sig_mix0_l.*CAL_l(:,1);
    results_SA_H0(pow_index,runindex) = yout_u(201:(200+Sam_Num))'*yout_u(201:(200+Sam_Num))/Sam_Num;
    results_DA_PED_H0(pow_index,runindex) = yout_u(201:(200+Sam_Num))'*yout_l(201:(200+Sam_Num))/Sam_Num;
    results_DA_OED_H0(pow_index,runindex) = norm((yout_u(201:(200+Sam_Num))+yout_l(201:(200+Sam_Num)))/2)^2/Sam_Num;
%     %% active cross-talk estimation & mitigation
%     Sam_end = 100+Sam_Num;
%     bhi = fir1(200,1/SampleNumberAve,'low');
%     y4_H1 = sig_mix1_u(1:Sam_end).*CAL_u(1:Sam_end,4);
%     temp_H1 = filter(bhi,1,y4_H1);
%     y_j4_hat_H1 = zeros(Sam_Num,1);
%     for jj=1:6
%         if jj==4
%             y_j4_hat_H1 = y_j4_hat_H1+temp_H1(101:Sam_end).*CAL_u(1:Sam_Num,4);
%         else
%         y_j4_hat_H1 = y_j4_hat_H1+real(temp_H1(101:Sam_end)...
%             .*exp(1j*(4-jj)/32*2*pi*(1:Sam_Num)')).*CAL_u(1:Sam_Num,jj);
%     %     y_j4_hat = y_j4_hat+real(sig_cache(1:5000,4)*PhaseShift(4)...
%     %         .*exp(1j*(4-jj)/32*2*pi*(1:5000)')).*CAL(1:5000,jj);
%         end
%     end
%     post_comp_H1 = (sig_mix1_u(1:Sam_Num)-y_j4_hat_H1);
%     results_SA_IC_H1(pow_index,runindex) = norm(filter(bhi,1,post_comp_H1.*CAL_u(1:Sam_Num,1)))^2/Sam_Num;
%         
%         
%     y4_H0 = sig_mix0_u(1:Sam_end).*CAL_u(1:Sam_end,4);
%     temp_H0 = filter(bhi,1,y4_H0);
%     y_j4_hat_H0 = zeros(Sam_Num,1);
%     for jj=1:6
%         if jj==4
%             y_j4_hat_H0 = y_j4_hat_H0+temp_H0(101:Sam_end).*CAL_u(1:Sam_Num,4);
%         else
%         y_j4_hat_H0 = y_j4_hat_H0+real(temp_H0(101:Sam_end)...
%             .*exp(1j*(4-jj)/32*2*pi*(1:Sam_Num)')).*CAL_u(1:Sam_Num,jj);
%     %     y_j4_hat = y_j4_hat+real(sig_cache(1:5000,4)*PhaseShift(4)...
%     %         .*exp(1j*(4-jj)/32*2*pi*(1:5000)')).*CAL(1:5000,jj);
%         end
%     end
%     post_comp_H0 = (sig_mix0_u(1:Sam_Num)-y_j4_hat_H0);
%     results_SA_IC_H0(pow_index,runindex) = norm(filter(bhi,1,post_comp_H0.*CAL_u(1:Sam_Num,1)))^2/Sam_Num;
%   
end

%
temp = sort(results_SA_H0(pow_index,:),'ascend');
TH_SA(pow_index) = temp(runtimes*0.9);
pd_SA(pow_index) = sum(results_SA_H1(pow_index,:)>TH_SA(pow_index))/runtimes;

temp = sort(results_DA_PED_H0(pow_index,:),'ascend');
TH_DA_PED(pow_index) = temp(runtimes*0.9);
pd_DA_PED(pow_index) = sum(results_DA_PED_H1(pow_index,:)>TH_DA_PED(pow_index))/runtimes;

temp = sort(results_DA_OED_H0(pow_index,:),'ascend');
TH_DA_OED(pow_index) = temp(runtimes*0.9);
pd_DA_OED(pow_index) = sum(results_DA_OED_H1(pow_index,:)>TH_DA_OED(pow_index))/runtimes;

% temp = sort(results_SA_IC_H0(pow_index,:),'ascend');
% TH_SA_IC(pow_index) = temp(runtimes*0.9);
% pd_SA_IC(pow_index) = sum(results_SA_IC_H1(pow_index,:)>TH_SA_IC(pow_index))/runtimes;
end
% %%
% NNN = 17e2:18e2;
% figure
% subplot(211)
% plot(sig_mix0_u(NNN));
% hold on
% plot(y_j4_hat_H0(NNN))
% hold on;
% ylim([-5,5])
% grid on
% legend('true','est.')
% subplot(212)
% plot(sig_mix0_u(NNN)-y_j4_hat_H0(NNN))
% ylim([-5,5])
% grid on
%%

% test1 = temp_H0(100+NNN).*CAL_u(NNN,4)...
%     +real(temp_H0(100+NNN).*exp(1j*(4-1)/32*2*pi*(NNN)')).*CAL_u(NNN,1);
% 
% test11 = temp_H0(100+NNN).*CAL_u(NNN,4);
% test12 = real(temp_H0(100+NNN).*exp(1j*(4-1)/32*2*pi*(NNN)'));
% 
% test2 = real(sqrt(2)*sig_cache(NNN,4)*PhaseShift(4)).*CAL_u(NNN,4)...
%     +real(sqrt(2)*sig_cache(NNN,4)*PhaseShift(4).*exp(1j*(4-1)/32*2*pi*(NNN)')).*CAL_u(NNN,1);
% 
% test21 = real(sqrt(2)*sig_cache(NNN,4)*PhaseShift(4)).*CAL_u(NNN,4);
% test22 = real(sqrt(2)*sig_cache(NNN,4)*PhaseShift(4).*exp(1j*(4-1)/32*2*pi*(NNN)'));
% 
% figure
% plot(test12);hold on
% plot(test22);hold on
% plot(y_j4_hat_H0(NNN));hold on
% plot(sig_mix0_u(NNN));
% legend('test','est','true')
% %%
% figure
% plot(temp_H0(100+NNN))
% hold on
% plot(sqrt(2)*real(sig_cache(NNN,4)*PhaseShift(4)));hold on
% % plot(y4_H0(NNN))
% legend('hat','true')
%%

figure(1)
plot(pow_pool,pd_SA,'o--','linewidth',3,'markersize',8);hold on
plot(pow_pool,pd_DA_PED,'o--','linewidth',3,'markersize',8);hold on
plot(pow_pool,pd_DA_OED,'o--','linewidth',3,'markersize',8);hold on
% plot(pow_pool,pd_SA_IC,'o--','linewidth',3,'markersize',8);hold on

grid on
legend('single-arm (ED0)','double-arm, proposed (ED1) ','double-arm, optimal (ED2)')
%% test distribution of ED0
Sam_Num = 5e3;
pow_index = 3;
s_pow = 10^(pow_pool(pow_index)/10)*(1+1/32*(6-1));
% s_pow = 10^(pow_pool(pow_index)/10)*0.5;
n_pow = 1/32*(6-1);
% n_pow = 0;
N_prime = Sam_Num/32;

% H0 emprical
[f,xi] = ksdensity(results_SA_H0(pow_index,:));
mu = mean(results_SA_H0(pow_index,:));
sigma = sqrt(var(results_SA_H0(pow_index,:)));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(99)
% plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
plot(xi,f,'--','linewidth',2);hold on

% H0 theoretical
mu = n_pow;
sigma = sqrt(2*n_pow^2/(1*N_prime));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(99)
plot(pdfxdata, normpdf(pdfxdata,mu,sigma),'linewidth',2);hold on
grid on

% H1 emprical
[f,xi] = ksdensity(results_SA_H1(pow_index,:));
mu = mean(results_SA_H1(pow_index,:));
sigma = sqrt(var(results_SA_H1(pow_index,:)));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(99)
% plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
plot(xi,f,'--','linewidth',2);hold on


% H1 theoretical
mu = s_pow+n_pow;
sigma = sqrt((2*(s_pow+n_pow)^2)/(1*N_prime));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(99)
plot(pdfxdata, normpdf(pdfxdata,mu,sigma),'linewidth',2);hold on
grid on
%% test distribution of ED1
pow_index = 3;
s_pow = 10^(pow_pool(pow_index)/10)*1;
n_pow = 1/32*(6-1);
N_prime = Sam_Num/32;

% H0 emprical
[f,xi] = ksdensity(results_DA_PED_H0(pow_index,:));
mu = mean(results_DA_PED_H0(pow_index,:));
sigma = sqrt(var(results_DA_PED_H0(pow_index,:)));
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
[f,xi] = ksdensity(results_DA_PED_H1(pow_index,:));
mu = mean(results_DA_PED_H1(pow_index,:));
sigma = sqrt(var(results_DA_PED_H1(pow_index,:)));
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
SINR_test = s_pow/n_pow;
TH_test = qfuncinv(0.1)*n_pow/sqrt(N_prime);
pd_test = qfunc((qfuncinv(0.1)-sqrt(N_prime)*SINR_test)/(sqrt((2*SINR_test^2+2*SINR_test+1))));
%
%% test distribution of ED2
pow_index = 3;
s_pow = 10^(pow_pool(pow_index)/10)*(1+0.5/32*(6-1));
% s_pow = 10^(pow_pool(pow_index)/10)*0.5;
n_pow = 0.5/32*(6-1);
% n_pow = 0;
N_prime = Sam_Num/32;

% H0 emprical
[f,xi] = ksdensity(results_DA_OED_H0(pow_index,:));
mu = mean(results_DA_OED_H0(pow_index,:));
sigma = sqrt(var(results_DA_OED_H0(pow_index,:)));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(99)
% plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
plot(xi,f,'--','linewidth',2);hold on

% H0 theoretical
mu = n_pow;
sigma = sqrt(2*n_pow^2/(1*N_prime));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(99)
plot(pdfxdata, normpdf(pdfxdata,mu,sigma),'linewidth',2);hold on
grid on

% H1 emprical
[f,xi] = ksdensity(results_DA_OED_H1(pow_index,:));
mu = mean(results_DA_OED_H1(pow_index,:));
sigma = sqrt(var(results_DA_OED_H1(pow_index,:)));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(99)
% plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
plot(xi,f,'--','linewidth',2);hold on


% H1 theoretical
mu = s_pow+n_pow;
sigma = sqrt((2*(s_pow+n_pow)^2)/(1*N_prime));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(99)
plot(pdfxdata, normpdf(pdfxdata,mu,sigma),'linewidth',2);hold on
grid on
%% Single-arm active cross-talk mitigation
pow_index = 3;
% H0 emprical
[f,xi] = ksdensity(results_SA_IC_H0(pow_index,:));
mu = mean(results_SA_IC_H0(pow_index,:));
sigma = sqrt(var(results_SA_IC_H0(pow_index,:)));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(98)
% plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
plot(xi,f,'--','linewidth',2);hold on

% H1 emprical
[f,xi] = ksdensity(results_SA_IC_H1(pow_index,:));
mu = mean(results_SA_IC_H1(pow_index,:));
sigma = sqrt(var(results_SA_IC_H1(pow_index,:)));
pdfxdata = linspace(mu-4*sigma,mu+4*sigma,100);
figure(98)
% plot(pdfxdata, normpdf(pdfxdata,mu,sigma));hold on
plot(xi,f,'--','linewidth',2);hold on
grid on
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
%% theoretical SINR verification
Sam_Num = Sam_Num;
SampleNumberAve_ED0 = 32;
SampleNumberAve_ED1 = 32;
SampleNumberAve_ED2 = 32;
N_prime_ED0 = Sam_Num/(SampleNumberAve_ED0*1.1);
N_prime_ED1 = Sam_Num/SampleNumberAve_ED1;
N_prime_ED2 = Sam_Num/(SampleNumberAve_ED0*1.1);

blk_num = 6;
pow_pool_fine = -10-linspace(0,20,100);
% pow_pool_fine = -18;
for pow_index = 1: 100
    s_pow_ED0(pow_index) = 10^(pow_pool_fine(pow_index)/10)*(1+1/SampleNumberAve_ED0*(blk_num-1));
    s_pow_ED1(pow_index) = 10^(pow_pool_fine(pow_index)/10)*1;
    s_pow_ED2(pow_index) = 10^(pow_pool_fine(pow_index)/10)*(1+0.5/SampleNumberAve_ED2*(blk_num-1));
    
    n_pow_ED0 = 1/SampleNumberAve_ED0*(blk_num-1)*1.1;
    n_pow_ED1 = 1/SampleNumberAve_ED1*(blk_num-1)*1.1;
    n_pow_ED2 = 0.5/SampleNumberAve_ED2*(blk_num-1)*1.1;
    
    SINR_ED0(pow_index) = s_pow_ED0(pow_index)/n_pow_ED0;
    SINR_ED1(pow_index) = s_pow_ED1(pow_index)/n_pow_ED1;
    SINR_ED2(pow_index) = s_pow_ED2(pow_index)/n_pow_ED2;
    
%     pd_theo_ED0(pow_index) = qfunc((qfuncinv(0.1)-sqrt(Sam_Num/2)*SINR(pow_index))/(1+SINR(pow_index)));
    
%     TH_theo_ED1(pow_index) = qfuncinv(0.1)*n_pow/sqrt(N_prime);
%     pd_theo_ED1(pow_index) = qfunc((-s_pow(pow_index)+TH_theo_ED1(pow_index))/(sqrt((2*s_pow(pow_index)^2+2*s_pow(pow_index)*n_pow+n_pow^2)/N_prime)));
    
%     pd_theo_ED2(pow_index) = qfunc((qfuncinv(0.1)-sqrt(N_prime/2)*SINR(pow_index))/(1+2*SINR(pow_index)));
    pd_theo_ED0_new(pow_index) = qfunc((qfuncinv(0.1)-sqrt(N_prime_ED0/2)*SINR_ED0(pow_index))/(1+SINR_ED0(pow_index)));
    pd_theo_ED1_new(pow_index) = qfunc((qfuncinv(0.1)-sqrt(N_prime_ED1)*SINR_ED1(pow_index))/(sqrt(1+2*SINR_ED1(pow_index)+2*SINR_ED1(pow_index)^2)));
    pd_theo_ED2_new(pow_index) = qfunc((qfuncinv(0.1)-sqrt(N_prime_ED2/2)*SINR_ED2(pow_index))/(1+SINR_ED2(pow_index)));

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

