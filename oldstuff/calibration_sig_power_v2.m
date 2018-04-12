%  05/24/2015
% initially we assume calibration signal has infinite strong power that
% second order terms of interferers are negligible. Here we tune such power
% to test how much performance will deteriorate.

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
M_num=12;
M_range=linspace(1,12,M_num);

LK_range= [8];

Pc_range=[1e0,1e-2,1e-4,1e-5,1e-6];
%Pc_range=1e6;
cal_pow_num=length(Pc_range);
pd=zeros(length(Pc_range),M_num);

record_table = zeros(cal_pow_num,M_num,runtimes_sim);

for runindex=1:runtimes_sim
    if mod(runindex,10)==0
        runindex/10
        CAL_perm=zeros(P,N);
        for mm=1:M
            CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
        end
    end
    start_point = randi(L-P-1,blk_num,1);
    sig_cache = zeros(P,blk_num);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii):start_point(ii)+P-1,randi(blk_num));
    end
    for cal_pow_index = 1:cal_pow_num
        Pc = Pc_range(cal_pow_index);
        sig_mix_perm_raw = zeros(P,1);
        for ii=1:blk_num
            sig_mix_perm_raw = sig_mix_perm_raw + 1/2/sqrt(Pc)*(sig_cache(:,ii).*sqrt(sig_pow(ii))).^2+...
            (sig_cache(:,ii)*sqrt(sig_pow(ii))).*CAL_perm(:,ii);
        end
        for mm=1:M_num
            M_sweep=M_range(mm);
            stage = 1;
            record_vec=zeros(1,8);
            sig_mix_perm = sig_mix_perm_raw(1:N*M_sweep);
            for ii=1:blk_num
                result_noproc_sim_perm(stage,ii) = tagging_v3(sig_mix_perm,CAL_perm(1:N*M_sweep,ii),N,M_sweep);
            end
            for stage=2:9
                [tempvalue,tempindex] = max(result_noproc_sim_perm(stage-1,:));
                record_vec(stage-1)=tempindex;
                sig_mix_perm = sig_mix_perm-1/2/sqrt(Pc)*(sig_cache(1:N*M_sweep,tempindex).*sqrt(sig_pow(tempindex))).^2-...
                    (sig_cache(1:N*M_sweep,tempindex)*sqrt(sig_pow(tempindex))).*CAL_perm(1:N*M_sweep,tempindex);
                for ii=1:blk_num
                    result_noproc_sim_perm(stage,ii) = tagging_v3(sig_mix_perm,CAL_perm(1:N*M_sweep,ii),N,M_sweep);
                end
            end
            record_table(cal_pow_index,mm,runindex)=sum((record_vec<=LK_range(1)))==LK_range(1);
        end
    end
end

%% detection performance
pd = zeros(cal_pow_num,M_num);
for cal_pow_index=1:cal_pow_num
    for mm=1:M_num
        pd(cal_pow_index,mm) = sum(squeeze(record_table(cal_pow_index,mm,:)))/runtimes_sim;
    end
end
%% figure
figure
plot_setting
plot(M_range*8,pd','-x');hold all
xlabel('Number of Frames')
ylabel('Probability of Detection')
legend('Pc = +\infty','Pc = -20dBm','Pc = -40dBm','Pc = -50dBm','Pc = -60dBm')
grid on
xlim([10,100])
ylim([0,1])
