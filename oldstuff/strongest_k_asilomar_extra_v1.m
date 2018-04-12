%  10/23/2015
%  an approximation has been used to theoretically calculate prob of
%  identification Kf strongest interferers.
%  such approximation take acoount error event of 
%  \hat{P}_{Kf+1} > \hat{P}_{Kf} and \hat{P}_{Kf+2} > \hat{P}_{Kf}

%  besides, here we sweep M. Theoretical approx. works bad with small M,
%  which means theoretical value can not predict poor performance regime.
%  it maybe because with small M, gaussian approx. does not match well
clear;clc;%clf;close all

blk_num = 31;
M = 1e2;
N = 31;
SampleNumberAve = 32;
stat_num = 2e2;
P = M*N;
L = M*SampleNumberAve*stat_num;

for ii = 1:blk_num
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(SampleNumberAve,1));
%     sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);
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
    %sig_i(:,ii) = zeros(size(sig_r(:,ii)));
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);

end
%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL_noperm=LowerRate_v2(cal0,P);
CAL_perm=zeros(P,N);
for mm=1:M
    CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
end
%% tagging
M_sweep_num=10;
M_sweep_range=linspace(5,50,M_sweep_num);
runtimes_sim = 1e4;
result_noproc_sim_perm=zeros(runtimes_sim,M,blk_num);
for runindex=1:runtimes_sim
    if mod(runindex,100)==1
        runindex/100
        sig_pow_dB = rand(blk_num,1)*60-60;
        sig_pow_dB = sort(sig_pow_dB,'descend');
        for kk=6:31
            if sig_pow_dB(kk)>(sig_pow_dB(5)-3)
                sig_pow_dB(kk) = -60;
            end
        end
        sig_pow = 10.^(sig_pow_dB/10);
        sig_pow = sig_pow./(sum(sig_pow));
        
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
index_cache = zeros(blk_num,runtimes_sim);
K_num=1;
K_range=[5];
k_largest_correct_flag=zeros(M_sweep_num,runtimes_sim,K_num);

for mm=1:M_sweep_num
    for runindex=1:runtimes_sim
        temp = squeeze(result_noproc_sim_perm(runindex,mm,:));
        [sortvalue,sortindex]= sort(temp,'descend');
        %index_cache(:,runindex) = sortindex;
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
%clf;close all
figure(1)
plot_setting
plot(M_sweep_range,k_largest_pd,'x');hold on
xlabel('M')
ylabel('Probability of Detection of k Largest Interferers')


