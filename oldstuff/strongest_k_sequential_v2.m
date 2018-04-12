%  05/25/2015
%  strongest K_T interferer detection performance. 
%  Both RRC and Rec pulas shape has been tested here.
clear;clc;clf;close all

load('cal_pow_test.mat');
%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL_noperm=LowerRate_v2(cal0,P);
CAL_perm=zeros(P,N);
for mm=1:M
    CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
end
%% tagging (sim.1 results)
M_num=10;
M_range=linspace(1,10,M_num);

LK_range= [6,7,8,9,10];
LK_num = length(LK_range);
pd=zeros(LK_num,M_num);
record_table=zeros(runtimes_sim,M_num,LK_num);
    
for runindex=1:runtimes_sim
    if mod(runindex,10)==0
        runindex/10
    end
    start_point = randi(L-P-1,blk_num,M);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii,mm):start_point(ii,mm)+P-1,randi(blk_num));
    end
    sig_mix_perm_raw = zeros(P,1);

    for ii=1:blk_num
        sig_mix_perm_raw = sig_mix_perm_raw + sig_cache(:,ii).*CAL_perm(:,ii)*sqrt(sig_pow(ii));
    end
    
    for mm=1:M_num
        M_sweep=M_range(mm);
        sig_mix_perm = sig_mix_perm_raw(1:N*M_sweep);
        record_vec = zeros(1,10);
        stage=1;
        for ii=1:blk_num
            result_noproc_sim_perm(stage,ii) = tagging_v3(sig_mix_perm,CAL_perm(1:N*M_sweep,ii),N,M_sweep);
        end
        for stage=2:11
            [tempvalue,tempindex] = max(result_noproc_sim_perm(stage-1,:));
            record_vec(stage-1)=tempindex;
            sig_mix_perm = sig_mix_perm-sig_cache(1:N*M_sweep,tempindex).*CAL_perm(1:N*M_sweep,tempindex)*sqrt(sig_pow(tempindex));
            for ii=1:blk_num
                result_noproc_sim_perm(stage,ii) = tagging_v3(sig_mix_perm,CAL_perm(1:N*M_sweep,ii),N,M_sweep);
            end
        end
        for kk=1:LK_num
            LK = LK_range(kk);
            record_table(runindex,mm,kk) = sum((record_vec<=LK))==LK;
        end
    end
end

%%
for mm=1:M_num
    for kk=1:LK_num
        pd(kk,mm) = sum(squeeze(record_table(:,mm,kk)))/runtimes_sim;
    end
end

%% figure
pd()
figure
plot_setting;
plot(M_range*10,pd','-x');hold all
xlabel('Number of Frames')
ylabel('Probability of Detection')
legend('Largest 6','Largest 7','Largest 8','Largest 9','Largest 10')
ylim([0,1])
