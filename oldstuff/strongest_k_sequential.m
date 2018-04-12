%  05/18/2015
% test for sequential detection.
% square waveform is considered here.
% this version uses real waveform.
clear;clc;clf;close all

load('matrix_iso_mean_63_rec.mat')
runtimes_sim = 1e3;

blk_num = 31;
M = 1e2;
N = 63;
SampleNumberAve = 40;
stat_num = 2e2;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = linspace(0,60,blk_num)'-90;
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);

for ii = 1:blk_num
    sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
    sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
%     upsam = SampleNumberAve;
%     symbols = fix(L*2/upsam);
%     
%     clearvars data temp_data
%     hmod = modem.pskmod('M', 4, 'InputType', 'integer');
%     hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
%     hpulse = design(hdesign);
%     data = randi(4,symbols,1)-1;
%     data = modulate(hmod, data);
%     data = upsample(data,upsam);
%     temp_data = conv(data,hpulse.Numerator);
%     unnormalized = real(temp_data(end-L+1:end));
%     sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     unnormalized = imag(temp_data(end-L+1:end));
%     sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     sig_i(:,ii) = zeros(size(sig_r(:,ii)));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))*sqrt(sig_pow(ii));

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
pd=zeros(3,M_num);
LK_range= [3,5,7];
for mm=1:M_num
    mm
    M_sweep=M_range(mm);
    record_table=zeros(runtimes_sim,10);
    result_noproc_sim_perm=zeros(10,blk_num);
    for runindex=1:runtimes_sim
        %CAL=CAL(:,randperm(N));
        start_point = randi(L-P-1,blk_num,M);
        for ii=1:blk_num
            sig_cache(:,ii) = sig(start_point(ii,mm):start_point(ii,mm)+P-1,ii);
        end
        sig_mix_perm = zeros(P,1);


        for ii=1:blk_num
            sig_mix_perm = sig_mix_perm + sig_cache(:,ii).*CAL_perm(:,ii);
        end
        stage=1;
        for ii=1:blk_num
            result_noproc_sim_perm(stage,ii) = tagging_v3(sig_mix_perm(1:N*M_sweep),CAL_perm(1:N*M_sweep,ii),N,M_sweep);
        end

        for stage=2:11
            [tempvalue,tempindex] = max(result_noproc_sim_perm(stage-1,:));
            record_table(runindex,stage-1)=tempindex;
            sig_mix_perm = sig_mix_perm-sig_cache(:,tempindex).*CAL_perm(:,tempindex);
            for ii=1:blk_num
                result_noproc_sim_perm(stage,ii) = tagging_v3(sig_mix_perm(1:N*M_sweep),CAL_perm(1:N*M_sweep,ii),N,M_sweep);
            end
        end
    end
    for kk=1:3
        LK=LK_range(kk);
        pd(kk,mm)=sum(sum(record_table<=LK,2)==LK)/runtimes_sim;
    end
end

%% figure
figure
plot_setting;
plot(M_range*10,pd','-x');hold all
xlabel('index')
ylabel('Probability of Detection')
legend('Largest 3','Largest 5','Largest 7')
ylim([0.3,1])
