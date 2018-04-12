%  05/25/2015
% strongest K detection
clear;clc;clf;close all

load('cal_pow_test_rec.mat');
%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL_noperm=LowerRate_v2(cal0,P);
CAL_perm=zeros(P,N);
for mm=1:M
    CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
end
%% tagging
M_sweep_num=10;
M_sweep_range=linspace(10,100,M_sweep_num);
runtimes_sim = 1e4;
result_noproc_sim_perm=zeros(runtimes_sim,M,blk_num);
for runindex=1:runtimes_sim
    if mod(runindex,100)==0
        runindex/100
        CAL_perm=zeros(P,N);
        for mm=1:M
            CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
        end
    end
    %CAL=CAL(:,randperm(N));
    start_point = randi(L-P-1,blk_num,M);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii,mm):start_point(ii,mm)+P-1,randi(blk_num));
    end
    sig_mix_perm = zeros(P,1);

    for ii=1:blk_num
        sig_mix_perm = sig_mix_perm + sig_cache(:,ii).*CAL_perm(:,ii)*sqrt(sig_pow(ii));
    end
    
    for ii=1:blk_num
        for mm=1:M_sweep_num
            M_sweep=M_sweep_range(mm);
            result_noproc_sim_perm(runindex,mm,ii) = tagging_v3(sig_mix_perm(1:N*M_sweep),CAL_perm(1:N*M_sweep,ii),N,M_sweep);
        end
    end
end

%% k largest detection performance
K_num=4;
K_range=[1,3,5,7];
k_largest_correct_flag=zeros(M_sweep_num,runtimes_sim,K_num);

for mm=1:M_sweep_num
    for runindex=1:runtimes_sim
        temp = squeeze(result_noproc_sim_perm(runindex,mm,:));
        [sortvalue,sortindex]= sort(temp,'descend');
        
        for kk=1:K_num
            K_sweep=K_range(kk);
            if sum(sort(sortindex(1:K_sweep))<=K_sweep)==K_sweep
                k_largest_correct_flag(mm,runindex,kk)=1;
            end
        end
    end
end
for mm=1:M_sweep_num
    for kk=1:K_num
        k_largest_pd(mm,kk)=sum(squeeze(k_largest_correct_flag(mm,:,kk)))/runtimes_sim;
    end
end
    
%% figure
figure
plot_setting
plot(M_sweep_range,k_largest_pd,'-x');
xlabel('M')
ylabel('Probability of Detection of k Largest Interferers')
legend('K=6','K=7','K=8','K=9','K = 10')

