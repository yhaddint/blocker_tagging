clear;clc;%clf;close all

blk_num = 31;
M = 1e2;
N = 31;
SampleNumberAve = 32;
stat_num = 2e2;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = linspace(0,60,blk_num)-60;
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);

for ii = 1:blk_num
    sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
    sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
    sig(:,ii) = (sig_r(:,ii)+1j*sig_i(:,ii))/sqrt(2);
%     upsam = SampleNumberAve;
%     symbols = fix(L*2/upsam);
%     
%     clearvars data temp_data
%     hmod = modem.qammod('M', 4, 'InputType', 'integer');
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
%     %sig_i(:,ii) = zeros(size(sig_r(:,ii)));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);

end

%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL_noperm=LowerRate_v2(cal0,P);
CAL_perm=zeros(P,N);
for mm=1:M
    CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
end
%% tagging
dft_mtx = dftmtx(N);

M_sweep = 100;
runtimes_sim = 1e3;
result_noproc_sim_perm=zeros(runtimes_sim,blk_num);
for runindex=1:runtimes_sim
    if mod(runindex,100)==1
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
        sig_mix_perm = sig_mix_perm + kron(ones(M,1),dft_mtx(:,ii)).*sig_cache(:,ii).*CAL_perm(:,ii)*sqrt(sig_pow(ii));
        %sig_mix_perm = sig_mix_perm + sig_cache(:,ii).*CAL_perm(:,ii)*sqrt(sig_pow(ii));
    end
    
    for ii=1:blk_num
        %decoding_seq = CAL_perm(1:N*M_sweep,ii);
        decoding_seq = kron(ones(M,1),conj(dft_mtx(:,3))).*CAL_perm(1:N*M_sweep,ii);
        result_noproc_sim_perm(runindex,ii) = tagging_v4(sig_mix_perm(1:N*M_sweep),decoding_seq,N,M_sweep);
    end
end
%%
for ii=1:blk_num
    pow_est(ii) = mean(squeeze(result_noproc_sim_perm(:,ii)));
end
10*log10(pow_est)