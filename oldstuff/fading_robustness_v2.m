%  12/18/2015
%  when signals have (block) fading, proposed method is supposed to be more robust
%  beacuse it takes more symbol duration to estimate power
%  
% Besides, the power of calibration signal itself is considered as finite
% value

clear;clc;%clf;close all

band_num = 8;
blk_num = 6;
M = 1e2;
N = 32;
SampleNumberAve = 128;

stat_num = 2e2;
P = M*N;
L = N*stat_num*M-P;

fading_flag = 0;
noise_power = 1e-3;
P_c = 10^(40/10);

rand('seed',4);
for ii = 1:blk_num
%     sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
%     sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
%     sig0(:,ii) = (sig_r(:,ii)+1j*sig_i(:,ii))/sqrt(2);
    upsam = SampleNumberAve;
    symbols = fix(L*2/upsam);
    
    clearvars data temp_data
    qam_level = 4;
    hmod = modem.qammod('M', qam_level, 'InputType', 'integer');
    hdesign  = fdesign.pulseshaping(upsam,'Square Root Raised Cosine');
    hpulse = design(hdesign);
    data = randi(qam_level,symbols,1)-1;
    data = modulate(hmod, data);
    data = upsample(data,upsam);
    temp_data = conv(data,hpulse.Numerator);
    sig0(:,ii) = temp_data./sqrt(sum(temp_data.*conj(temp_data))/length(temp_data));
end

%% Formulating multi-band signals
band_realization_num = 300;
F_pool = [3];
F_num = length(F_pool);
ave_num_pool = 1:1:24;
pd = zeros(length(ave_num_pool),F_num,band_realization_num);

P_ur_dB = -8;
P_ur=10^(P_ur_dB/10);

for band_realization_index=1:band_realization_num
    disp(band_realization_index/band_realization_num*100)
    % signal table (band | subband | power)
    sig_table = zeros(blk_num,3);
    % sig in which band
    sig_table(:,1) = randi(band_num,blk_num,1);
    % sig in which subband
    sig_table(:,2) = randi(N/2-1,blk_num,1);
    % power of sig
%     sig_pow_dB = [-3,-6,-18,-30,-30,-30];
%     sig_pow = 10.^(sig_pow_dB/10);
%     sig_table(1:6,3) = sig_pow';
    sig_pow_dB = [0.42075;0.42075];
%     sig_pow_dB = sort(rand(2,1)*8-8,'descend');
    sig_pow0 = 10.^(sig_pow_dB/10);
    sig_pow = sig_pow0./sum(sig_pow0)*(1-P_ur);
    sig_table(1:3,3) = [sig_pow;P_ur];

    for ff=1:max(F_pool)
        best_bin_ref(ff) = (sig_table(ff,1)-1)*N/2+sig_table(ff,2);
    end

    Nn = N;
    runtimes = 20;

    % fading_coeff = raylrnd(1,blk_num,runtimes);
    pd_counter = zeros(runtimes,length(F_pool),length(ave_num_pool));

    for runindex = 1:runtimes
        LL = N*max(ave_num_pool);
        y = zeros(LL,blk_num);
        CAL = randi(2,LL,band_num)*2-3;
        for ii=1:blk_num
            start_P = randi(L);
    %         sig_cache0 = sig(start_P:start_P+N-1,ii)*sqrt(10^(-randn));
            sig_cache0 = sig0(start_P:start_P+LL-1,ii);
            sig_cache1 = sig_cache0.*transpose(exp(1j*(rand*2*pi+(sig_table(ii,2)-1)*2*pi*(1:LL)/N)))*sqrt(sig_table(ii,3));
            if fading_flag
%                 chan = rayleighchan(1/6e6,100);
%                 sig_cache = filter(chan,sig_cache1);
                sig_cache = randn*real(sig_cache1)+1j*randn*imag(sig_cache1);
            else
                sig_cache = sig_cache1;
            end
            y(:,ii) = 2*real(sig_cache).*CAL(:,sig_table(ii,1))+abs(sig_cache).^2/sqrt(P_c);
        end
        sig_sum = sum(y,2);%+sqrt(noise_power)*(rand(N,1)*2-1);
        
        % fft-based periodogram estimation
        bins_bands_temp = zeros(max(ave_num_pool),N/2,band_num);
        for kk=1:band_num
            sig_sum_corr = sig_sum.*CAL(:,kk);
            for aa=1:max(ave_num_pool)
                periodo_gram = (abs(fft(sig_sum_corr((aa-1)*N+1:aa*N))/N).^2)';
                bins_bands_temp(aa,:,kk) = periodo_gram(1:N/2);
            end
        end
        
        for ave_index = 1:length(ave_num_pool)
            ave_num = ave_num_pool(ave_index);
            for kk=1:band_num
                bins_bands(:,kk) = sum(squeeze(bins_bands_temp(1:ave_num,:,kk)),1)'/ave_num;
            end

            % detection
            
            [sortval,sortindex] = sort((bins_bands(:)),'descend');
            for ff=1:length(F_pool)
                F=F_pool(ff);
%                 bestbin(runindex,:)=sortindex(1:F);
                if sort(sortindex(1:F))'==sort(best_bin_ref(1:F))
                    pd_counter(runindex,ff,ave_index) = pd_counter(runindex,ff,ave_index)+1;
                end
            end
        end
    end

    for ave_index = 1:length(ave_num_pool)
        for ff=1:length(F_pool)
            pd(ave_index,ff,band_realization_index) = sum(pd_counter(:,ff,ave_index))/runtimes;
        end
    end
end
%% average over different realization of frequency settings
for ave_index=1:length(ave_num_pool)
    for ff=1:F_num
        pd_diff_freq(ave_index,ff) = mean(pd(ave_index,ff,:));
    end
end


F_num = length(F_pool);
target_pool = 0.5:0.05:0.9;
ydata = zeros(length(target_pool),F_num);
for ff=1:F_num
    for tt=1:length(target_pool)
        target = target_pool(tt);
        if min(pd_diff_freq(:,ff))>target
            ydata(tt,ff)=1;
        elseif max(pd_diff_freq(:,ff))<target
            ydata(tt,ff)=0;
        else
            xx = min(find(pd_diff_freq(:,ff)>target));
            ydata(tt,ff)= xx-1+1/(pd_diff_freq(xx,ff)-pd_diff_freq(xx-1,ff))*(target-pd_diff_freq(xx-1,ff));
        end
    end
end
figure
plot(target_pool,ydata*N,'-o');hold on

%%
for ff=2
    target_pool = 0.5:0.025:0.9;
    pd_cand = linspace(min(pd_diff_freq(:,ff)),max(pd_diff_freq(:,ff)),1e3);
    sample_int = interp1(pd_diff_freq(:,ff),ave_num_pool*N,pd_cand);
    for ii=1:length(target_pool)
        [bestscore_index,best_index] = min(abs(pd_cand-target_pool(ii)));
        ydata(ii,ff) = sample_int(best_index);
    end
end
figure
plot(target_pool,ydata,'-o');hold on

%% identification figure
% figure(4)
% plot_setting();
% plot(ave_num_pool*Nn,pd_diff_freq,'x--');
% hold on
% grid on
% xlabel('Number of Samples')
% ylabel('Ident. Prob.')
% ylim([0,1])
% % legend('F = 1', 'F = 2', 'F = 3')
%%
% figure
% plot(interp1(ave_num_pool,pd_diff_freq,linspace(2,32,20),'spline'),linspace(2,32,20)*N,'-o');
% hold on
%%
% figure
% plot(linspace(min(pd_diff_freq),max(pd_diff_freq),20),interp1(pd_diff_freq,ave_num_pool*N,linspace(min(pd_diff_freq),max(pd_diff_freq),20),'spline'));
% hold on
%%
figure
plot(ave_num_pool*Nn,pd_diff_freq);
%%
% figure(5)
% plot_setting();
% plot(pd_diff_freq,ave_num_pool*Nn);
% hold on
% grid on
% ylabel('Required Samples')
% xlabel('Targeted Ident. Prob.')
% xlim([0.5,1])
% % legend('F = 1', 'F = 2', 'F = 3')
%% spectrum figure
% figure(1)
% for kk=1:band_num
%     temp = squeeze(bins_bands(:,:,kk));
%     ydata = mean(abs(temp),1);
%     subplot(band_num/2,2,kk)
%     plot(10*log10(ydata));
%     grid on
% end

% % figure
% figure(3)
% 
% temp = squeeze(y_bins_bands(:,:,1));
% ydata = mean(abs(temp),1);
% plot(10*log10(ydata));
% grid on
