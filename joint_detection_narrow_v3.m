%  04/14/2015
% test for joint estimation.
% RRC waveform is considered here.(square waveform is in version 2)
% this version uses real waveform.
clear;clc;clf;close all
load('matrix_iso_mean_63_v3.mat')

runtimes_sim = 2e3;
runtimes_analy = 1e4;

blk_num = 30;
matrix_proc_sep = inv(matrix_iso_mean(1:blk_num,1:blk_num));
M = 1e2;
N = 63;
SampleNumberAve = 20;
stat_num = 3e3;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = linspace(0,60,blk_num)'-90;
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);

for ii = 1:blk_num
    upsam = SampleNumberAve;
    symbols = fix(L*2/upsam);
    
    clearvars data temp_data
    hmod = modem.pskmod('M', 4, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(4,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
    unnormalized = real(temp_data(end-L+1:end));
    sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
    unnormalized = imag(temp_data(end-L+1:end));
    sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
    sig_i(:,ii) = zeros(size(sig_r(:,ii)));
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))*sqrt(sig_pow(ii));

end

%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL=LowerRate_v2(cal0,P);

% CAL=zeros(P,N);
% for mm=1:M
%     CAL((mm-1)*N+1:mm*N,:)=CAL_temp((mm-1)*N+1:mm*N,randperm(N));
% end
%% process matrix
% N = 63, R = 20
% mus=0.3039;
% mun=0.0113;
% sigmas=0.2256;
% sigman=0.235e-3;

% N = 31, R = 20
mus=0.6257;
mun=0.01279;
sigmas=0.6716;
sigman=0.5e-3;
diag_comp=1/(mus-mun);
offdiag_comp=-mun/(mus-mun)/(mus+(blk_num-1)*mun);
matrix_proc=ones(ii)*offdiag_comp+diag(ones(1,ii)*(diag_comp));

%% tagging (sim.1 results)
for runindex=1:runtimes_sim
    %CAL=CAL(:,randperm(N));
    start_point = randi(L-P-1,blk_num,M);
    for ii=1:blk_num
        for mm=1:M
            sig_cache1((mm-1)*N+1:mm*N,ii) = sig(start_point(ii,mm):start_point(ii,mm)+N-1,ii);
        end
        sig_cache2(:,ii) = sig(start_point(ii,mm):start_point(ii,mm)+P-1,ii);
    end
    sig_mix1 = zeros(P,1);
    sig_mix2 = zeros(P,1);

    for ii=1:blk_num
        sig_mix1 = sig_mix1 + sig_cache1(:,ii).*CAL(:,ii);
        sig_mix2 = sig_mix2 + sig_cache2(:,ii).*CAL(:,ii);
    end
    
    for ii=1:blk_num
        result_noproc_sim1(runindex,ii) = tagging_v3(sig_mix1,CAL(:,ii),N,M)/mus;
        result_noproc_sim2(runindex,ii) = tagging_v3(sig_mix2,CAL(:,ii),N,M)/mus;

    end
    
    result_proc_sim1(runindex,:) = (matrix_proc*result_noproc_sim1(runindex,:)')'*mus;
    result_proc_sim2(runindex,:) = (matrix_proc*result_noproc_sim2(runindex,:)')'*mus;

    result_proc_sep(runindex,:) = (matrix_proc_sep*result_noproc_sim1(runindex,:)')'*mus;
end

%% analytical results

mu_matrix=ones(blk_num,blk_num)*mun;
mu_matrix=mu_matrix+diag(ones(1,blk_num)*(mus-mun));

sig_get = mu_matrix*sig_pow;

result_noproc_analy=zeros(runtimes_analy,blk_num);
sig_var=sig_pow.*(matrix_iso_var(1:blk_num,1:blk_num)*sig_pow)/M;

for ii=1:blk_num
%     sigma_matrix=diag(ones(1,blk_num))*sigman;
%     sigma_matrix(ii,ii)=sigmas;
%     for jj=1:blk_num
%         for kk=1:blk_num
%             if (jj~=kk) 
%                if jj==ii || kk==ii
%                    sigma_matrix(jj,kk)=2*mus*mun;
%                else
%                    sigma_matrix(jj,kk)=2*mun*mun;
%                end
%             end
%         end
%     end
    result_noproc_analy(:,ii)=(sig_get(ii)+randn(runtimes_analy,1)*sqrt(sig_var(ii)))/mus;
end

for rr=1:runtimes_analy
    result_proc_analy(rr,:) = (matrix_proc*result_noproc_analy(rr,:)'*mus)';
end
%% detection
% TH=-45;
% pd_sim1 = sum(result_proc_sim1 > (10^(TH/10)))/runtimes_sim;
% pd_sim2 = sum(result_proc_sim2 > (10^(TH/10)))/runtimes_sim;
% 
% pd_analy = sum(result_proc_analy > (10^(TH/10)))/runtimes_analy;
% 
% % plot
% figure(1)
% plot_setting
% plot(sig_pow_dB,pd_analy,'r');hold on
% plot(sig_pow_dB,pd_sim1,'rx');hold on
% plot(sig_pow_dB,pd_sim2,'ro');hold on
% 
% xlabel('Interferer Power (dBm)')
% ylabel('Prob being Detected')
% legend('Analy.','Sim.1','Sim.2')
% %
% figure(2)
% plot(sort(result_proc_analy(1:runtimes_sim,15)),'b');hold on
% plot(sort(result_proc_sim1(1:runtimes_sim,15)),'r')
% plot(sort(result_proc_sim2(1:runtimes_sim,15)),'c')
%
figure(3)
plot(10*log10(var(result_noproc_analy(:,:))),'b');hold on
plot(10*log10(var(result_proc_analy(:,:))),'bx');hold on

plot(10*log10(var(result_noproc_sim1(:,:))),'r');hold on
plot(10*log10(var(result_proc_sim1(:,:))),'rx');hold on

plot(10*log10(var(result_proc_sep(:,:))),'g');hold on
title('Var')

figure(4)
clearvars temp
temp = mean(result_noproc_analy(:,:));
plot(10*log10(temp.*(temp>0)),'b');hold on
temp = mean(result_proc_analy(:,:));
plot(10*log10(temp.*(temp>0)),'bx');hold on

temp = mean(result_noproc_sim1(:,:));
plot(10*log10(temp.*(temp>0)),'r');hold on
temp = mean(result_proc_sim1(:,:));
plot(10*log10(temp.*(temp>0)),'rx');hold on

temp = mean(result_proc_sep(:,:));
plot(10*log10(temp.*(temp>0)),'go');hold on

% temp = mean(result_noproc_sim2(:,:));
% plot(10*log10(temp.*(temp>0)),'r--');hold on
% temp = mean(result_proc_sim2(:,:));
% plot(10*log10(temp.*(temp>0)),'ro');hold on

plot(sig_pow_dB,'k--')
ylim([-90,-30]) 
title('Mean')