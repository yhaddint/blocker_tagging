%  04/14/2015
% test for sequential detection.
% square waveform is considered here.
% this version uses real waveform.
clear;clc;clf;close all

load('matrix_iso_mean_63_rec.mat')
runtimes_sim = 1e4;
runtimes_analy = 1e4;

blk_num = 31;
M = 1e2;
N = 63;
SampleNumberAve = 40;
stat_num = 2e3;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = linspace(0,60,blk_num)'-90;
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);

for ii = 1:blk_num
    sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
    sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
    sig_i(:,ii) = zeros(size(sig_r(:,ii)));
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))*sqrt(sig_pow(ii));
end


%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL_noperm=LowerRate_v2(cal0,P);

%% process matrix
% mus=0.284;
% mun=0.0115;
% sigmas=0.1854;
% sigman=0.5e-3;

mus=0.5006;
mun=0.008178;
sigmas=0.1830;
sigman=0.33e-3;
diag_comp=1/(mus-mun);
offdiag_comp=-mun/(mus-mun)/(mus+(blk_num-1)*mun);
matrix_proc_perm=ones(ii)*offdiag_comp+diag(ones(1,ii)*(diag_comp));
matrix_proc_noperm = inv(matrix_iso_mean(1:blk_num,1:blk_num));
%% tagging (sim.1 results)
for runindex=1:runtimes_sim
    % permutation of code
    CAL_perm=zeros(P,N);
    for mm=1:M
        CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
    end

    %CAL=CAL(:,randperm(N));
    start_point = randi(L-P-1,blk_num,M);
    for ii=1:blk_num
        for mm=1:M
            sig_cache1((mm-1)*N+1:mm*N,ii) = sig(start_point(ii,mm):start_point(ii,mm)+N-1,ii);
        end
        sig_cache2(:,ii) = sig(start_point(ii,mm):start_point(ii,mm)+P-1,ii);
    end
    sig_mix1_perm = zeros(P,1);
    sig_mix2_perm = zeros(P,1);
    sig_mix1_noperm = zeros(P,1);
    sig_mix2_noperm = zeros(P,1);

    for ii=1:blk_num
        sig_mix1_perm = sig_mix1_perm + sig_cache1(:,ii).*CAL_perm(:,ii);
        sig_mix1_noperm = sig_mix1_noperm + sig_cache1(:,ii).*CAL_noperm(:,ii);
        sig_mix2_perm = sig_mix2_perm + sig_cache2(:,ii).*CAL_perm(:,ii);
        sig_mix2_noperm = sig_mix2_noperm + sig_cache2(:,ii).*CAL_noperm(:,ii);

    end
    
    for ii=1:blk_num
        result_noproc_sim1_perm(runindex,ii) = tagging_v3(sig_mix1_perm,CAL_perm(:,ii),N,M);
        result_noproc_sim1_noperm(runindex,ii) = tagging_v3(sig_mix1_noperm,CAL_noperm(:,ii),N,M);
        result_noproc_sim2_perm(runindex,ii) = tagging_v3(sig_mix2_perm,CAL_perm(:,ii),N,M);
        result_noproc_sim2_noperm(runindex,ii) = tagging_v3(sig_mix2_noperm,CAL_noperm(:,ii),N,M);

    end
    result_proc_sim1_perm(runindex,:) = (matrix_proc_perm*result_noproc_sim1_perm(runindex,:)')';
    result_proc_sim1_noperm(runindex,:) = (matrix_proc_noperm*result_noproc_sim1_noperm(runindex,:)')';
    result_proc_sim2_perm(runindex,:) = (matrix_proc_perm*result_noproc_sim2_perm(runindex,:)')';
    result_proc_sim2_noperm(runindex,:) = (matrix_proc_noperm*result_noproc_sim1_noperm(runindex,:)')';

end

%% analytical results

% mu_matrix=ones(blk_num,blk_num)*mun;
% mu_matrix=mu_matrix+diag(ones(1,blk_num)*(mus-mun));
% 
% sig_get = mu_matrix*sig_pow;
% 
% result_noproc_analy=zeros(runtimes_analy,blk_num);
% for ii=1:blk_num
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
%     sig_var(ii)=sig_pow'*sigma_matrix*sig_pow/M;
%     result_noproc_analy(:,ii)=(sig_get(ii)+randn(runtimes_analy,1)*sqrt(sig_var(ii)))/mus;
% end
% 
% for rr=1:runtimes_analy
%     result_proc_analy(rr,:) = (matrix_proc_perm*result_noproc_analy(rr,:)'*mus)';
% end
%% detection
% TH=-50;
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
% legend('Analy.','Sim.')
% %
% figure(2)
% subplot(211)
% plot(sort(result_proc_analy(1:runtimes_sim,15)),'b');hold on
% plot(sort(result_proc_sim1(1:runtimes_sim,15)),'r')
% plot(sort(result_proc_sim2(1:runtimes_sim,15)),'c')
% subplot(212)
% plot(sort(result_noproc_analy(1:runtimes_sim,15)),'b');hold on
% plot(sort(result_noproc_sim1_perm(1:runtimes_sim,15)),'r')
% plot(sort(result_noproc_sim2_perm(1:runtimes_sim,15)),'c')

%%
% figure(3)
% plot(10*log10(var(result_noproc_analy(:,:))),'b');hold on
% plot(10*log10(var(result_proc_analy(:,:))),'b--');hold on
% 
% plot(10*log10(var(result_noproc_sim1_perm(:,:))),'r');hold on
% plot(10*log10(var(result_proc_sim1(:,:))),'rx');hold on
% 
% plot(10*log10(var(result_noproc_sim2_perm(:,:))),'c');hold on
% plot(10*log10(var(result_proc_sim2(:,:))),'cx');hold on
% title('Var')
%%
figure(4)
plot_setting;
clearvars temp
% temp = mean(result_noproc_analy(:,:));
% plot(10*log10(temp.*(temp>0)),'b');hold on
% temp = mean(result_proc_analy(:,:));
% plot(10*log10(temp.*(temp>0)),'bx');hold on

temp = mean(result_noproc_sim1_perm(:,:));
plot(10*log10(temp.*(temp>0)),'bs');hold on
temp = mean(result_proc_sim1_perm(:,:));
plot(10*log10(temp.*(temp>0)),'bo');hold on

% temp = mean(result_noproc_sim1_noperm(:,:));
% plot(10*log10(temp.*(temp>0)),'r');hold on
% temp = mean(result_proc_sim1_noperm(:,:));
% plot(10*log10(temp.*(temp>0)),'rx');hold on

temp = mean(result_noproc_sim2_perm(:,:));
plot(10*log10(temp.*(temp>0)),'gs');hold on
temp = mean(result_proc_sim2_perm(:,:));
plot(10*log10(temp.*(temp>0)),'go');hold on

% temp = mean(result_noproc_sim2_noperm(:,:));
% plot(10*log10(temp.*(temp>0)),'m--');hold on
% temp = mean(result_proc_sim2_noperm(:,:));
% plot(10*log10(temp.*(temp>0)),'mo');hold on
plot(sig_pow_dB,'kx')

legend('Without Biasing Compensation (Random Frame)','With Biasing Compensation (Random Frame)','Without Biasing Compensation (Consecutive Frame)','With Biasing Compensation (Consecutive Frame)','True Power Value')

ylim([-90,-30])
xlim([1,31])
xlabel('Interferer Index')
ylabel('Power (dBm)')
title('Mean')
%%
% figure(5)
% 
% % temp = mean(result_noproc_analy(:,:));
% % plot(10*log10(temp.*(temp>0)),'b');hold on
% % temp = mean(result_proc_analy(:,:));
% % plot(10*log10(temp.*(temp>0)),'bx');hold on
% 
% plot(mean(result_noproc_sim1_perm(:,:)),'c');hold on
% plot(mean(result_proc_sim1_perm(:,:)),'cx');hold on
% 
% plot(mean(result_noproc_sim1_noperm(:,:)),'r');hold on
% plot(mean(result_proc_sim1_noperm(:,:)),'rx');hold on
% 
% 
% plot(mean(result_noproc_sim2_perm(:,:)),'g--');hold on
% plot(mean(result_proc_sim2_perm(:,:)),'go');hold on
% 
% plot(mean(result_noproc_sim2_noperm(:,:)),'m--');hold on
% plot(mean(result_proc_sim2_noperm(:,:)),'mo');hold on
% 
% plot(sig_pow,'k--')
% %ylim([-90,-30])
% title('Mean')