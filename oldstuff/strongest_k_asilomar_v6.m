%  10/29/2015
% strongest K detection. sweeping N
clear;clc;%clf;close all
fig_data=zeros(4,1,4);
N_range=[31,63,127,255];
for nn=1:4
    
N=N_range(nn);    
blk_num = 31;
M = 50;
SampleNumberAve = 32;
stat_num = 2e2;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = linspace(0,90,blk_num)'-90;
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);

for ii = 1:blk_num
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
result_noproc_sim_perm=zeros(runtimes_sim,blk_num);
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
        result_noproc_sim_perm(runindex,ii) = tagging_v3(sig_mix_perm(1:N*M),CAL_perm(1:N*M,ii),N,M);
    end
end

%% k largest detection performance
K_num=4;
K_range=[1,3,5,7];
k_largest_correct_flag=zeros(runtimes_sim,K_num);


for runindex=1:runtimes_sim
    temp = result_noproc_sim_perm(runindex,:);
    [sortvalue,sortindex]= sort(temp,'descend');

    for kk=1:K_num
        K_sweep=K_range(kk);
        if sum(sort(sortindex(1:K_sweep))<=K_sweep)==K_sweep
            k_largest_correct_flag(runindex,kk)=1;
        end
    end
end

for kk=1:K_num
    k_largest_pd(kk)=sum(k_largest_correct_flag(:,kk))/runtimes_sim;
    fig_data(nn,kk) = k_largest_pd(kk);
end

    

end
%% figure
figure
plot_setting
plot(N_range,transpose(squeeze(fig_data(:,:,1))),'b^');hold on
plot(N_range,transpose(squeeze(fig_data(:,:,2))),'gx');hold on
plot(N_range,transpose(squeeze(fig_data(:,:,3))),'ro');hold on
plot(N_range,transpose(squeeze(fig_data(:,:,4))),'cs');hold on

% plot(N_range,(squeeze(fig_data(:,:,1))),'b-x');hold on
% plot(N_range,(squeeze(fig_data(:,:,2))),'b-x');hold on
% plot(N_range,(squeeze(fig_data(:,:,3))),'b-x');hold on
% plot(N_range,(squeeze(fig_data(:,:,4))),'b-x');hold on

%plot(N_range,squeeze(fig_data(:,:,3)),'-x');hold on
%plot(N_range,squeeze(fig_data(:,:,4)),'-x');hold on
xlabel('N')
ylabel('Identification Probability')
legend('F=1','F=3','F=5','F=7')
%%
blk_num = 31;
sig_pow_dB = linspace(0,60,blk_num)'-60;
sig_pow_dB = sort(sig_pow_dB,'descend');
sig_pow = 10.^(sig_pow_dB/10);
K_num=4;
K_range=[1,3,5,7];
N_range=[31,63,127,255];
upsilon_pool = [0.7899,0.4332,0.2416,0.1234];
mu_pool = [0.007818,0.008773,0.006034,0.003453];
M_sweep_range = 50;
ident_prob_approx = zeros(length(K_range),length(N_range));

for kk=1:K_num
    kf = K_range(kk);
    for nn=1:length(N_range)
        upsilon = upsilon_pool(nn);
        mu = mu_pool(nn);
        EP1 = sig_pow(kf)*upsilon+(sum(sig_pow)-sig_pow(kf))*mu;
        EP2 = sig_pow(kf+1)*upsilon+(sum(sig_pow)-sig_pow(kf+1))*mu;
        sigmaP1 = sqrt(EP1^2*2/M_sweep_range);
        sigmaP2 = sqrt(EP2^2*2/M_sweep_range);
        ident_prob_approx(kk,nn) = 1-qfunc((EP1-EP2)/sqrt(sigmaP1^2+sigmaP2^2));%-integral(fun22,0,5);
    end
end
% figure
figure(1)
plot(N_range,ident_prob_approx');hold on
legend('F=1 (Sim.)','F=3 (Sim.)','F=5 (Sim.)','F=7 (Sim.)','F=1 (Aprx.)','F=3 (Aprx.)','F=5 (Aprx.)','F=7 (Aprx.)')
ylim([0,1])

