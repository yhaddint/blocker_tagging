%  06/19/2015
% case study. its simumation setting details is listted in the slides
% here power of interferers are determinstic
clear;clc;clf;close all


runtimes_sim = 1e3;

blk_num = 31;
M = 30;
N = 31;
SampleNumberAve = 32;
stat_num = 5e2;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = linspace(-3,-93,blk_num)';
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
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);

end

%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL_noperm=LowerRate_v2(cal0,P);
CAL_perm=zeros(P,N);
for mm=1:M
    CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
end
%% tagging: first stage
P_TH1 = 10^(-1.4);
P_TH2 = 10^(-2);

result_stage1 = zeros(runtimes_sim,blk_num);
result_stage2 = zeros(runtimes_sim,blk_num);

for runindex=1:runtimes_sim
    start_point = randi(L-P-1,blk_num,M);
    for ii=1:blk_num
        sig_cache(:,ii) = sig(start_point(ii,mm):start_point(ii,mm)+P-1,ii).*sqrt(sig_pow(ii));
    end
    sig_mix_perm = zeros(P,1);

    for ii=1:blk_num
        sig_mix_perm = sig_mix_perm + sig_cache(:,ii).*CAL_perm(:,ii);
    end

    for ii=1:blk_num
        result_stage1(runindex,ii) = tagging_v3(sig_mix_perm,CAL_perm(1:N*M,ii),N,M);
    end
    
    tag_table_s1(runindex,:) = result_stage1(runindex,:)>=P_TH1;
    
    % stage 2
    sig_mix_perm = zeros(P,1);
    for ii=1:blk_num
        if tag_table_s1(runindex,ii)~=1
            sig_mix_perm = sig_mix_perm + sig_cache(:,ii).*CAL_perm(:,ii);
        end
    end

    for ii=1:blk_num
        result_stage2(runindex,ii) = tagging_v3(sig_mix_perm,CAL_perm(1:N*M,ii),N,M);
    end
    tag_table_s2(runindex,:) = result_stage2(runindex,:)>=P_TH2;
end
%% detection
pdpfa_s1 = mean(tag_table_s1(:,:));
pdpfa_s2 = mean(tag_table_s2(:,:));
pdpfa_all = mean(tag_table_s1(:,:)+tag_table_s2(:,:));
%% figure
figure
plot_setting;
plot(pdpfa_s1,'b');hold all
plot(pdpfa_s2,'r');hold all
plot(pdpfa_all,'c');hold all
xlabel('index')
ylabel('Probability of Detection')

