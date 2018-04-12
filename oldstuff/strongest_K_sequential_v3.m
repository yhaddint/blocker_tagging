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
runtimes_sim = 5e3;
LK_range= [6,8,10];
%LK_range = 6;
LK_num = length(LK_range);

record_table=zeros(runtimes_sim,20,LK_num);
    
for runindex=1:runtimes_sim
    if mod(runindex,20)==0
        runindex/20
        CAL_perm=zeros(P,N);
        for mm=1:M
            CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
        end
    end
    start_point = randi(L-P-1,blk_num,M);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii,mm):start_point(ii,mm)+P-1,randi(blk_num));
    end
    sig_mix_perm_raw = zeros(P,1);

    for ii=1:blk_num
        sig_mix_perm_raw = sig_mix_perm_raw + sig_cache(:,ii).*CAL_perm(:,ii)*sqrt(sig_pow(ii));
    end
    for kk=1:LK_num
        LK=LK_range(kk);
        M_range=LK:LK:100;
        M_num=length(M_range);
        for mm=1:M_num
            M_sweep=mm;
            sig_mix_perm = sig_mix_perm_raw(1:N*M_sweep);
            record_vec = zeros(1,LK);
            stage=1;
            for ii=1:blk_num
                result_noproc_sim_perm(stage,ii) = tagging_v3(sig_mix_perm,CAL_perm(1:N*M_sweep,ii),N,M_sweep);
            end
            for stage=2:LK+1
                [tempvalue,tempindex] = max(result_noproc_sim_perm(stage-1,:));
                record_vec(stage-1)=tempindex;
                sig_mix_perm = sig_mix_perm-sig_cache(1:N*M_sweep,tempindex).*CAL_perm(1:N*M_sweep,tempindex)*sqrt(sig_pow(tempindex));
                for ii=1:blk_num
                    result_noproc_sim_perm(stage,ii) = tagging_v3(sig_mix_perm,CAL_perm(1:N*M_sweep,ii),N,M_sweep);
                end
            end
            record_table(runindex,mm,kk) = sum((record_vec<=LK))==LK;
        end
    end
end

%%
figure
plot_setting

for kk=1:LK_num
    LK=LK_range(kk);
    M_range=LK:LK:100;
    M_num=length(M_range);
    pd=zeros(1,M_num);
    for mm=1:M_num
        pd(mm) = sum(squeeze(record_table(:,mm,kk)))/runtimes_sim;
    end
    plot(M_range,pd);hold on
end
xlabel('Number of Frames')
ylabel('Probability of Detection')
legend('Largest 6','Largest 8','Largest 10')
