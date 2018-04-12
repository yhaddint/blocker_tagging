%  05/23/2015
%  using inv(U) to compensate biasing. I'd like to show biasing in mean can
%  be entirely compensate while variance is still very large.
%  Both RRC and Rec will be tested while I believe only one will be
%  presented in the paper to make it clear.
clear;clc;clf;close all

runtimes_sim = 1e4;

blk_num = 31;
M = 1e2;
N = 63;
SampleNumberAve = 40;
stat_num = 2e3;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = linspace(0,40,blk_num)'-70;
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

%% process matrix
mus=0.5006;
mun=0.008178;
% sigmas=0.1830;
% sigman=0.33e-3;
diag_comp=1/(mus-mun);
offdiag_comp=-mun/(mus-mun)/(mus+(blk_num-1)*mun);
matrix_proc_perm=ones(ii)*offdiag_comp+diag(ones(1,ii)*(diag_comp));

mus=0.5006;
mun=1/N;
% sigmas=0.1830;
% sigman=0.33e-3;
diag_comp=1/(mus-mun);
offdiag_comp=-mun/(mus-mun)/(mus+(blk_num-1)*mun);
matrix_proc_rc = ones(ii)*offdiag_comp+diag(ones(1,ii)*(diag_comp));
%% tagging (sim.1 results)
for runindex=1:runtimes_sim
    % permutation of code
    CAL_perm=zeros(P,N);
    for mm=1:M
        CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
    end
    CAL_rc=randi(2,M*N,blk_num)*2-3;
    
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
    sig_mix1_rc = zeros(P,1);
    sig_mix2_rc = zeros(P,1);

    for ii=1:blk_num
%         sig_mix1_perm = sig_mix1_perm + sig_cache1(:,ii).*CAL_perm(:,ii);
        %sig_mix1_rc = sig_mix1_rc + sig_cache1(:,ii).*CAL_rc(:,ii);
        sig_mix2_perm = sig_mix2_perm + sig_cache2(:,ii).*CAL_perm(:,ii);
        sig_mix2_rc = sig_mix2_rc + sig_cache2(:,ii).*CAL_rc(:,ii);

    end
    
    for ii=1:blk_num
%         result_noproc_sim1_perm(runindex,ii) = tagging_v3(sig_mix1_perm,CAL_perm(:,ii),N,M);
        %result_noproc_sim1_rc(runindex,ii) = tagging_v3(sig_mix1_rc,CAL_rc(:,ii),N,M);
        result_noproc_sim2_perm(runindex,ii) = tagging_v3(sig_mix2_perm,CAL_perm(:,ii),N,M);
        result_noproc_sim2_rc(runindex,ii) = tagging_v3(sig_mix2_rc,CAL_rc(:,ii),N,M);

    end
%     result_proc_sim1_perm(runindex,:) = (matrix_proc_perm*result_noproc_sim1_perm(runindex,:)')';
    %result_proc_sim1_rc(runindex,:) = (matrix_proc_perm*result_noproc_sim1_rc(runindex,:)')';
    result_proc_sim2_perm(runindex,:) = (matrix_proc_perm*result_noproc_sim2_perm(runindex,:)')';
    result_proc_sim2_rc(runindex,:) = (matrix_proc_rc*result_noproc_sim2_rc(runindex,:)')';

end

%% mean test with pn codes
figure
plot_setting
clearvars temp

% temp = mean(result_noproc_sim1_perm(:,:));
% plot(10*log10(temp.*(temp>0)),'bs');hold on
% temp = mean(result_proc_sim1_perm(:,:));
% plot(10*log10(temp.*(temp>0)),'bo');hold on

temp = mean(result_noproc_sim2_perm(:,:));
plot(10*log10(temp.*(temp>0)),'bs');hold on
temp = mean(result_proc_sim2_perm(:,:));
plot(10*log10(temp.*(temp>0)),'bo');hold on

plot(sig_pow_dB,'kx')

legend('Without Biasing Compensation (Random Frame)','With Biasing Compensation (Random Frame)','Without Biasing Compensation (Consecutive Frame)','With Biasing Compensation (Consecutive Frame)','True Power Value')

ylim([-75,-30])
xlim([1,31])
xlabel('Interferer Index')
ylabel('Power (dBm)')
title('Mean')
% mean test with random codes

clearvars temp

% temp = mean(result_noproc_sim1_rc(:,:));
% plot(10*log10(temp.*(temp>0)),'bs');hold on
% temp = mean(result_proc_sim1_rc(:,:));
% plot(10*log10(temp.*(temp>0)),'bo');hold on
% 
temp = mean(result_noproc_sim2_rc(:,:));
plot(10*log10(temp.*(temp>0)),'r^');hold on
temp = mean(result_proc_sim2_rc(:,:));
plot(10*log10(temp.*(temp>0)),'rp');hold on

plot(sig_pow_dB,'kx')

legend('Without Biasing Compensation (Random Frame)','With Biasing Compensation (Random Frame)','Without Biasing Compensation (Consecutive Frame)','With Biasing Compensation (Consecutive Frame)','True Power Value')

ylim([-70,-30])
xlim([1,31])
xlabel('Interferer Index')
ylabel('Power (dBm)')
title('Mean')
%% variance test PN codes
figure
plot_setting
% temp = sqrt(var(result_noproc_sim1_perm(:,:)));
% plot(10*log10(temp.*(temp>0)),'r');hold on
% temp = sqrt(var(result_proc_sim1_perm(:,:)));
% plot(10*log10(temp.*(temp>0)),'rx');hold on
% % 

temp = sqrt(var(result_noproc_sim2_perm(:,:)));
plot(10*log10(temp.*(temp>0)),'bs');hold on
temp = sqrt(var(result_proc_sim2_perm(:,:)));
plot(10*log10(temp.*(temp>0)),'bo');hold on
% variance test random codes

% temp = sqrt(var(result_noproc_sim1_perm(:,:)));
% plot(10*log10(temp.*(temp>0)),'r');hold on
% temp = sqrt(var(result_proc_sim1_perm(:,:)));
% plot(10*log10(temp.*(temp>0)),'rx');hold on
% % 

temp = sqrt(var(result_noproc_sim2_rc(:,:)));
plot(10*log10(temp.*(temp>0)),'r^');hold on
temp = sqrt(var(result_proc_sim2_rc(:,:)));
plot(10*log10(temp.*(temp>0)),'rp');hold on
