%  05/24/2015
% initially we assume calibration signal has infinite strong power that
% second order terms of interferers are negligible. Here we tune such power
% to test how much performance will deteriorate.

%  it doesn't work this time!
%  fix it tomorrow!!!!!
clear;clc;clf;close all


runtimes_sim = 1e3;

blk_num = 31;
M = 1e2;
N = 63;
SampleNumberAve = 40;
stat_num = 200;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = linspace(0,60,blk_num)'-90;
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);

for ii = 1:blk_num
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig_i(:,ii) = zeros(size(sig_r(:,ii)));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))*sqrt(sig_pow(ii));
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
    unnormalized = real(temp_data(end-L+1:end));
    sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
    unnormalized = imag(temp_data(end-L+1:end));
    sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
    sig_i(:,ii) = zeros(size(sig_r(:,ii)));
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
end

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

LK_range= [7];

%Pc_range=[1e2,1e1,1e0,1e-1,1e-2,1e-3];
Pc_range=1e6;
cal_pow_num=length(Pc_range);
pd=zeros(length(Pc_range),M_num);

record_table = zeros(cal_pow_num,M_num,runtimes_sim);

for runindex=1:runtimes_sim
    if mod(runindex,10)==0
        runindex/10
    end
    start_point = randi(L-P-1,blk_num,1);
    sig_cache = zeros(P,blk_num);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii):start_point(ii)+P-1,randi(blk_num));
    end
    for cal_pow_index = 1:cal_pow_num
        Pc = Pc_range(cal_pow_index);
        sig_mix_perm = zeros(P,1);
        for ii=1:blk_num
            sig_mix_perm = sig_mix_perm + 1/2/sqrt(Pc)*(sig_cache(:,ii).*sqrt(sig_pow(ii))).^2+...
            (sig_cache(:,ii)*sqrt(sig_pow(ii))).*CAL_perm(:,ii);
        end
        for mm=1:M_num
            M_sweep=M_range(mm);
            stage = 1;
            record_vec=zeros(1,10);
            for ii=1:blk_num
                result_noproc_sim_perm(stage,ii) = tagging_v3(sig_mix_perm(1:N*M_sweep),CAL_perm(1:N*M_sweep,ii),N,M_sweep);
            end
            for stage=2:11
                [tempvalue,tempindex] = max(result_noproc_sim_perm(stage-1,:));
                record_vec(stage-1)=tempindex;
                sig_mix_perm = sig_mix_perm-(1/2/sqrt(Pc)*(sig_cache(:,tempindex).*sqrt(sig_pow(tempindex))).^2+...
                    (sig_cache(:,tempindex)*sqrt(sig_pow(tempindex))).*CAL_perm(:,tempindex));
                for ii=1:blk_num
                    result_noproc_sim_perm(stage,ii) = tagging_v3(sig_mix_perm(1:N*M_sweep),CAL_perm(1:N*M_sweep,ii),N,M_sweep);
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
plot(M_range*10,pd','-x');hold all
xlabel('index')
ylabel('Probability of Detection')
legend('Pc=1e2','Pc=1e1','Pc=1e0','Pc=1e-1','Pc=1e-2','Pc=1e-3')

