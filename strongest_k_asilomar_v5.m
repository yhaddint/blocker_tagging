%  09/02/2015
% strongest K detection. sweeping symbol duration period
clear;clc;clf;close all
R_range=[8,16,32,64,128];
fig_data=zeros(4,4);
for R_index=1:4
clearvars -except R_range R_index fig_data
blk_num = 31;
M = 1e2;
N = 63;
SampleNumberAve = R_range(R_index);
stat_num = 2e2;
P = M*N;
L = M*SampleNumberAve*stat_num;

sig_pow_dB = linspace(0,90,blk_num)'-90;
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
    sig(:,ii) = (sig_r(:,ii)+sig_i(:,ii))/sqrt(2);

end
%% calibration sequences
cal0=PNgenerator_v5(N,N,1);
CAL_noperm=LowerRate_v2(cal0,P);
CAL_perm=zeros(P/2,N);
for mm=1:M/2
    CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
end
%% tagging
M_sweep_num=1;
M_sweep_range=50;
runtimes_sim = 5e3;
result_noproc_sim_perm=zeros(runtimes_sim,M,blk_num);
for runindex=1:runtimes_sim
    if mod(runindex,100)==0
        runindex/100
        CAL_perm=zeros(P/2,N);
        for mm=1:M/2
            CAL_perm((mm-1)*N+1:mm*N,:)=CAL_noperm((mm-1)*N+1:mm*N,randperm(N));
        end
    end
    %CAL=CAL(:,randperm(N));
    start_point = randi(L-P-1,blk_num,M);
    for ii=1:blk_num
        for mm=1:M_sweep_range
        sig_cache((mm-1)*N+1:mm*N,ii) = sig(start_point(ii,mm):start_point(ii,mm)+N-1,randi(blk_num));
        end
    end
    sig_mix_perm = zeros(P/2,1);

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
    end
end
    

fig_data(:,R_index) =  k_largest_pd;
end
%% figure
clf;close all
figure(1)
plot(R_range(1:4),fig_data','x');hold on
grid on

% theoretical approximation
%upsilon_pool = [0.4378,0.7165,0.9084,0.9754];
%mu_pool = [0.01918,0.01019,0.003997,0.001834];

% upsilon_pool = [0.235,0.4319,0.7102,0.9056];
% mu_pool = [0.01239,0.009269,0.004853,0.00175];
% 
% ident_prob_approx = zeros(length(K_range),length(R_range));
% for kk=1:K_num
%     kf = K_range(kk);
%     for rr=1:length(R_range)
%         
%         upsilon = upsilon_pool(rr);
%         mu = mu_pool(rr);
%         EP1 = sig_pow(kf)*upsilon+(sum(sig_pow)-sig_pow(kf))*mu;
%         EP2 = sig_pow(kf+1)*upsilon+(sum(sig_pow)-sig_pow(kf+1))*mu;
%         EP3 = sig_pow(kf+2)*upsilon+(sum(sig_pow)-sig_pow(kf+2))*mu;
%         sigmaP1 = sqrt(EP1^2*2/M_sweep_range);
%         sigmaP2 = sqrt(EP2^2*2/M_sweep_range);
%         sigmaP3 = sqrt(EP3^2*2/M_sweep_range);
% 
%         fun1 = @(x) exp(-0.5*((x-EP1).^2/sigmaP1^2+(x-EP2).^2/sigmaP2^2))/sqrt(2*pi*sigmaP1^2);
%         fun21 = @(x) qfunc((x-EP2)/sigmaP2).*exp(-0.5*(x-EP1).^2/sigmaP1^2)./sqrt(2*pi*sigmaP1^2);
%         fun22 = @(x) qfunc((x-EP3)/sigmaP3).*exp(-0.5*(x-EP1).^2/sigmaP1^2)./sqrt(2*pi*sigmaP1^2);
%         fun31 = @(x) (1/4*exp(-((x-EP2)/sigmaP2).^2/2)+1/4*exp(-((x-EP2)/sigmaP2).^2/2)).*exp(-0.5*(x-EP1).^2/sigmaP1^2)./sqrt(2*pi*sigmaP1^2);
% 
%         ident_prob_approx(kk,rr) = 1-integral(fun31,0,5);%-integral(fun22,0,5);
%     end
% end
% % figure
% figure(1)
% plot(R_range,ident_prob_approx');hold on
% legend('K=1 (Sim.)','K=3 (Sim.)','K=5 (Sim.)','K=7 (Sim.)','K=1 (Anal.)','K=3 (Anal.)','K=5 (Anal.)','K=7 (Anal.)')
% ylim([0,1])
%%
upsilon_pool = [0.235,0.4319,0.7102,0.9056];
mu_pool = [0.01239,0.009269,0.004853,0.00175];

ident_prob_approx = zeros(length(K_range),length(R_range));
for kk=1:K_num
    kf = K_range(kk);
    for rr=1:length(R_range)
        upsilon = upsilon_pool(rr);
        mu = mu_pool(rr);
        EP1 = sig_pow(kf)*upsilon+(sum(sig_pow)-sig_pow(kf))*mu;
        EP2 = sig_pow(kf+1)*upsilon+(sum(sig_pow)-sig_pow(kf+1))*mu;
        sigmaP1 = sqrt(EP1^2*2/M_sweep_range);
        sigmaP2 = sqrt(EP2^2*2/M_sweep_range);

        ident_prob_approx(kk,rr) = 1-qfunc((EP1-EP2)/sqrt(sigmaP1^2+sigmaP2^2));%-integral(fun22,0,5);
    end
end
% figure
figure(1)
plot(R_range,ident_prob_approx');hold on
legend('K=1 (Sim.)','K=3 (Sim.)','K=5 (Sim.)','K=7 (Sim.)','K=1 (Anal.)','K=3 (Anal.)','K=5 (Anal.)','K=7 (Anal.)')
ylim([0,1])