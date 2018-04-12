%  09/02/2015
%  precision of C2I2 estimator is simulated here.
%  We uses LFSR
%  Sim and Anal are compared.
%  In this version 1 we use rectangular pulse shape
clear;clc;clf;close all

runtimes_sim = 1e4;

blk_num = 31;
M = 30;
N = 31;
SampleNumberAve = 40;
stat_num = 2e2;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = linspace(0,30,blk_num)-30';
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);

for ii = 1:blk_num
    sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
    sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
end


%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL_noperm=LowerRate_v2(cal0,P);


%% tagging (sim.1 results)
for runindex=1:runtimes_sim
    % permutation of code
    CAL_perm=zeros(P,N);
    for mm=1:M
        CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
    end
    %CAL_rc=randi(2,M*N,blk_num)*2-3;
    
    start_point = randi(L-P-1,blk_num,M);
    bi=randi(blk_num,blk_num,1);
    for ii=1:blk_num
%         for mm=1:M
%             sig_cache1((mm-1)*N+1:mm*N,ii) = sig(start_point(ii,mm):start_point(ii,mm)+N-1,bi(ii))*sqrt(sig_pow(ii));
%         end
        sig_cache2(:,ii) = sig(start_point(ii,mm):start_point(ii,mm)+P-1,bi(ii))*sqrt(sig_pow(ii));
    end
    sig_mix1_perm = zeros(P,1);
    sig_mix2_perm = zeros(P,1);
    %sig_mix1_rc = zeros(P,1);
    %sig_mix2_rc = zeros(P,1);

    for ii=1:blk_num
%         sig_mix1_perm = sig_mix1_perm + sig_cache1(:,ii).*CAL_perm(:,ii);
        %sig_mix1_rc = sig_mix1_rc + sig_cache1(:,ii).*CAL_rc(:,ii);
        sig_mix2_perm = sig_mix2_perm + sig_cache2(:,ii).*CAL_perm(:,ii);
        %sig_mix2_rc = sig_mix2_rc + sig_cache2(:,ii).*CAL_rc(:,ii);

    end
    
    for ii=1:blk_num
%         result_noproc_sim1_perm(runindex,ii) = tagging_v3(sig_mix1_perm,CAL_perm(:,ii),N,M);
        %result_noproc_sim1_rc(runindex,ii) = tagging_v3(sig_mix1_rc,CAL_rc(:,ii),N,M);
        result_noproc_sim2_perm(runindex,ii) = tagging_v3(sig_mix2_perm,CAL_perm(:,ii),N,M);
        %result_noproc_sim2_rc(runindex,ii) = tagging_v3(sig_mix2_rc,CAL_rc(:,ii),N,M);

    end

end

%% mean test with LFSR codes
figure(99)
plot_setting
clearvars temp


temp = mean(result_noproc_sim2_perm(:,:));
plot(10*log10(temp.*(temp>0)),'bs');hold on

plot(sig_pow_dB,'kx')


ylim([-45,0])
xlim([1,31])
xlabel('Interferer Power dBm')
ylabel('Mean of Estimator (dBm)')
title('Mean')

%% variance test with LFSR codes
figure(98)
plot_setting
clearvars temp


temp = sqrt(var(result_noproc_sim2_perm(:,:)));
plot(10*log10(temp.*(temp>0)),'bs');hold on

plot(sig_pow_dB,'kx')


ylim([-45,0])
xlim([1,31])
xlabel('Interferer Power (dBm)')
ylabel('Standard Variance of Estimator (dBm)')
title('Var')
% mean test with random codes



%% anal mean
blk_num = 31;
sig_pow_dB = linspace(0,30,blk_num)'-30;
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);

mu_data = [0.7419,0.871;0.5006,0.7376];
upsilon_data = [0.009365,0.005203;0.008178,0.004418];
pow_est=zeros(31,2,2);
for nn=1:1
    for rr=1:1
        mu=mu_data(nn,rr);
        upsilon=upsilon_data(nn,rr);
        for ii=1:blk_num
            pow_est(ii,nn,rr) = sig_pow(ii)*mu+(sum(sig_pow)-sig_pow(ii))*upsilon;
        end
    end
end
figure(99)
for nn=1:1
    for rr=1:1
        mean_anal = squeeze(pow_est(:,nn,rr));
        plot(10*log10(squeeze(pow_est(:,nn,rr))));hold on
    end
end
%% anal var
var_anal = sqrt(2/M).*mean_anal;
figure(98)
plot(10*log10(var_anal));hold on
legend('','','Approx.')