%  09/02/2015
% strongest K detection. sweeping N
clear;clc;%clf;close all
fig_data=zeros(4,1,4);
N_range=[31,63,127,255];
for nn=1:4
    
N=N_range(nn);    
blk_num = 31;
M = 1e2;
SampleNumberAve = 40;
stat_num = 2e2;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = linspace(0,60,blk_num)'-90;
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);

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
M_sweep_num=1;
M_sweep_range=50;
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
    sig_cache = zeros(P,blk_num);
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
        fig_data(nn,mm,kk) = k_largest_pd(mm,kk);
    end
end
    

end
%% figure
figure
plot_setting
plot(N_range,transpose(squeeze(fig_data(:,:,1))),'b-x');hold on
plot(N_range,transpose(squeeze(fig_data(:,:,2))),'r-x');hold on
plot(N_range,transpose(squeeze(fig_data(:,:,3))),'g-x');hold on
plot(N_range,transpose(squeeze(fig_data(:,:,4))),'c-x');hold on

% plot(N_range,(squeeze(fig_data(:,:,1))),'b-x');hold on
% plot(N_range,(squeeze(fig_data(:,:,2))),'b-x');hold on
% plot(N_range,(squeeze(fig_data(:,:,3))),'b-x');hold on
% plot(N_range,(squeeze(fig_data(:,:,4))),'b-x');hold on

%plot(N_range,squeeze(fig_data(:,:,3)),'-x');hold on
%plot(N_range,squeeze(fig_data(:,:,4)),'-x');hold on
xlabel('N')
ylabel('Probability')
legend('K_F=1','K_F=3','K_F=5','K_F=7')

