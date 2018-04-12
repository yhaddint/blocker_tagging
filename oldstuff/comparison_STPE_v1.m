%  12/16/2015
%  test how many samples are required for STPE method in finding strongest
%  blocker.
%  




clear;clc;%clf;close all

band_num = 2;
blk_num = 16;
M = 1e2;
N = 32;
SampleNumberAve = 32;

stat_num = 2e2;
P = M*N;
L = N*stat_num;

noise_power = 1e-3;

rand('seed',1);
for ii = 1:blk_num
    sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
    sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
    sig0(:,ii) = (sig_r(:,ii)+1j*sig_i(:,ii))/sqrt(2);
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
%     sig(:,ii) = temp_data./sqrt(sum(temp_data.*conj(temp_data))/length(temp_data));

end

%% Formulating multi-band signals
band_realization_num = 1;
F_pool = [1 2 3];
F_num = length(F_pool);
ave_num_pool = [4];
pd = zeros(length(ave_num_pool),F_num,band_realization_num);

for band_realization_index=1:band_realization_num
%% Formulating multi-band signals
% signal table (band | subband | power)
sig_table = zeros(blk_num,3);
% sig in which band
sig_table(:,1) = randi(band_num,blk_num,1);
% sig in which subband
sig_table(:,2) = randi(N/2-1,blk_num,1);
% power of sig
sig_pow_dB = [-3,-6,-9,-20,-20];
sig_pow = 10.^(sig_pow_dB/10);
sig_table(1:5,3) = sig_pow';

for ii = 1:blk_num
    freq_mask = transpose(exp(1j*(sig_table(ii,2)-1)*2*pi*(1:length(sig0(:,ii)))/N));
    sig(:,ii) = sig0(:,ii).*freq_mask*sqrt(sig_table(ii,3));
end

runtimes = 2000;
%sig_conv = zeros(runtimes,L);
bins_bands = zeros(runtimes,N/2,band_num);
for runindex = 1:runtimes
    for kk=1:band_num
        y = zeros(P,blk_num);
        start_point = randi(P*M,1,blk_num);
        for ii=1:blk_num
            if sig_table(ii,1)==kk
            y(:,ii) = sig(start_point(ii):start_point(ii)+P-1,ii);
            end
        end
        sig_sum = sum(y,2)+sqrt(noise_power)*(rand(P,1)*2-1);
        periodo_gram = (fft(sig_sum(1:N))/N)';
        periodo_square = abs(periodo_gram).^2;
        bins_bands(runindex,:,kk) = [periodo_square(1),periodo_square(2:N/2)];
    end
end

%% ident. prob. calculation
best_bin_ref = zeros(1,max(F_pool));
for F_index = 1:F_num
    F = F_pool(F_index);
    for ff=1:F
        best_bin_ref(ff) = (sig_table(ff,1)-1)*N/2+sig_table(ff,2);
    end
    bestbin = zeros(runtimes,F);
    A = zeros(N/2,band_num);
    for ave_index = 1:length(ave_num_pool)
        clearvars bestscore bestbin
        pd_counter = 0;
        ave_num = ave_num_pool(ave_index);
        if ave_num == 1
            for ii=1:runtimes
               A = abs(squeeze(bins_bands(ii,:,:)));
               [sortval,sortindex] = sort((A(:)),'descend');
               bestbin(ii,:)=sortindex(1:F);
               if sort(sortindex(1:F))'==sort(best_bin_ref(1:F))
                    pd_counter = pd_counter+1;
               end
            end
        else
            for ii = 1:runtimes/ave_num
                for kk=1:band_num
                    A(:,kk) = mean(squeeze(bins_bands((ii-1)*ave_num+1:ii*ave_num,:,kk)));
                end
                [sortval,sortindex] = sort((A(:)),'descend');
                bestbin(ii,:)=sortindex(1:F);
                if sort(sortindex(1:F))'==sort(best_bin_ref(1:F))
                    pd_counter = pd_counter+1;
                end
            end
        end
        pd(ave_index,F_index,band_realization_index) = pd_counter/(runtimes/ave_num);
    end
end
end
%% average over different realization of frequency settings
for ave_index=1:length(ave_num_pool)
    for F_index = 1:F_num
        pd_diff_freq(ave_index,F_index) = mean(pd(ave_index,F_index,:));
    end
end
%% identification figure
figure
plot_setting();
plot(ave_num_pool*N*band_num,pd_diff_freq,'o-');
hold on
grid on
xlabel('Number of Samples')
ylabel('Ident. Prob.')
ylim([0,1])
legend('F = 1', 'F = 2', 'F = 3')
%% spectrum figure
figure
for kk=1:band_num
    temp = squeeze(bins_bands(:,:,kk));
    ydata = mean(abs(temp),1);
    subplot(band_num/2,2,kk)
    plot(10*log10(ydata));
    grid on
end

