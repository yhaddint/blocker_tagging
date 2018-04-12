% 05/27/2015
% estimation uncertainties are mutually uncorrealted?
%  Here is a simple test
clear;clc;clf;close all

runtimes_sim = 1e4;

blk_num = 31;
M = 1e2;
N = 63;
SampleNumberAve = 40;
stat_num = 500;
P = M*N;
L = M*SampleNumberAve*stat_num;

K = 10;
sig_pow = zeros(blk_num,K);
for kk=1:10
    sig_pow(1:kk,kk) = 1/kk;
end
sig_pow_dB = ones(blk_num,K)*(-60);
sig_pow = sig_pow + 10.^(sig_pow_dB/10);

for ii = 1:blk_num
    sig_r(:,ii) = kron(randi(4,M*stat_num,1)*2-5,ones(SampleNumberAve,1));
    sig_i(:,ii) = kron(randi(4,M*stat_num,1)*2-5,ones(SampleNumberAve,1));
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
    
%     upsam = SampleNumberAve;
%     symbols = fix(L*2/upsam);
%     
%     clearvars data temp_data
%     hmod = modem.qammod('M', 16, 'InputType', 'integer');
%     hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
%     hpulse = design(hdesign);
%     data = randi(16,symbols,1)-1;
%     data = modulate(hmod, data);
%     data = upsample(data,upsam);
%     temp_data = conv(data,hpulse.Numerator);
%     unnormalized = real(temp_data(end-L+1:end));
%     sig_r(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     unnormalized = imag(temp_data(end-L+1:end));
%     sig_i(:,ii) = unnormalized./sqrt(sum(unnormalized.^2)/L);
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);

end

%% tagging (sim.1 results)
result_noproc_sim_perm=zeros(runtimes_sim,blk_num,K);
for runindex=1:runtimes_sim
    if mod(runindex,100)==0
        runindex/100
    end
    %CAL=CAL(:,randperm(N));
    CAL = randi(2,P,blk_num)*2-3;
    start_point = randi(L-P-1,blk_num,M);
    sig_cache = zeros(P,blk_num);
    for ii=1:blk_num
        for mm=1:M
            sig_cache((mm-1)*N+1:mm*N,ii) = sig(start_point(ii,mm):start_point(ii,mm)+N-1,randi(blk_num));
        end
    end
    sig_mix_perm = zeros(P,K);

    for kk=1:10
        for ii=1:blk_num
            sig_mix_perm(:,kk) = sig_mix_perm(:,kk) + sig_cache(:,randi(blk_num)).*CAL(:,ii)*sqrt(sig_pow(ii,kk));
        end

        for ii=1:blk_num
            M_sweep=40;
            result_noproc_sim_perm(runindex,ii,kk) = tagging_v3(sig_mix_perm(1:N*M_sweep,kk),CAL(1:N*M_sweep,ii),N,M_sweep);
        end
    end
end

%%
ave = zeros(500,10);
count = 1;
for n1=11:31
    for n2=11:31
        if (n1~=n2)
            for kk=1:10
                E1(kk) = mean(squeeze(result_noproc_sim_perm(:,n1,kk)));
                v1(kk)= var(squeeze(result_noproc_sim_perm(:,n1,kk)));
                E2(kk) = mean(squeeze(result_noproc_sim_perm(:,n2,kk)));
                v2(kk)= var(squeeze(result_noproc_sim_perm(:,n2,kk)));
                E12(kk) = mean(squeeze(result_noproc_sim_perm(:,n1,kk)).*squeeze(result_noproc_sim_perm(:,n2,kk)));
            end
        ave(count,:) = (E12-E1.*E2)./(sqrt(v1.*v2));
        count=count+1;
        end
    end
end
ydata = sum(ave)/(21*20);
%%
plot(ydata)
%plot(1:K,ave/(20)/(19));