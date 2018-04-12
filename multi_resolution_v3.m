%  12/01/2015
%  test multi-resolution idea
%  Pseudo-spectrum is used. 
%  force R(0) to be zero to improve identification prob.
%  avoid using abs() to cancel out imagine part of bins
%  

clear;clc;clf;close all

blk_num = 31;
M = 1e2;
N = 32;
SampleNumberAve = 32;

stat_num = 2e2;
P = M*N;
L = N*stat_num;

sig_pow_dB = [-30,0,-30];
sig_pow = 10.^(sig_pow_dB/10);
rand('seed',1);
for ii = 1:blk_num
    sig_r(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
    sig_i(:,ii) = kron(randi(2,M*stat_num,1)*2-3,ones(N,1));
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
%     sig(:,ii) = temp_data./sqrt(sum(temp_data.*conj(temp_data))/length(temp_data));

end
freq_mask1 = transpose(exp(1j*6*2*pi*(1:L*M)/N));
sig(:,2) = sig(:,2).*freq_mask1;

freq_mask2 = transpose(exp(1j*2*2*pi*(1:L*M)/N));
sig(:,3) = sig(:,3).*freq_mask2;
%% a simple test with 3 interferers.
%  First two are in the same band but different subband, the third is in
%  another band

Nn = 4*N;
L = N;
runtimes = 5e3;
%sig_conv = zeros(runtimes,L);
bins_bands = zeros(runtimes,L/2,2);
for runindex = 1:runtimes
    phi = exp(1j*2*pi*rand(1,3));
    %phi = ones(1,3);
    sig_sum = zeros(P,1);
    CAL = randi(2,P,3)*2-3;
    sig_cache = zeros(P,3);
    start_point = randi(P*M,1,3);
    for ii=1:3
        sig_cache(:,ii) = sig(start_point(ii):start_point(ii)+P-1,ii);
        if ii ~= 3
            y(:,ii) = 2*real(sig_cache(:,ii).*phi(ii)).*CAL(:,1)*sqrt(sig_pow(ii));
        else
            y(:,ii) = 2*real(sig_cache(:,ii).*phi(ii)).*CAL(:,2)*sqrt(sig_pow(ii));
        end
        sig_sum = sum(y,2);
    end
%     
    for ii=1:2
        sig_sum_corr = sig_sum.*CAL(:,ii);
%         sig_conv(runindex,1)=0;
%         for ll=2:L
%             sig_conv(runindex,ll) = sum(sig_sum_corr(1:Nn).*(sig_sum_corr((1:Nn)+ll-1)))/Nn;
%         end
%         bins_bands_uncom = fft(sig_conv(runindex,:))/L;
        sig_conv(1)=0;
        for ll=2:L
            sig_conv(ll) = sum(sig_sum_corr(1:Nn).*(sig_sum_corr((1:Nn)+ll-1)))/Nn;
        end
        bins_bands_uncom = fft(sig_conv)/L;
        
        bins_bands(runindex,:,ii) = ([bins_bands_uncom(1),bins_bands_uncom(2:L/2)+fliplr(bins_bands_uncom(L/2+2:L))]);
        %bins_bands(runindex,:,ii) = ([bins_bands_uncom(1),2*bins_bands_uncom(2:L/2)]);
    end
    
    % watch for interference of sig1 to signal 3
%     y_corr = y(:,3).*CAL(:,2);
%     y_conv(1) = 0;
%     for ll=2:L
%         y_conv(ll) = sum(y_corr(1:Nn).*(y_corr((1:Nn)+ll-1)))/Nn;
%     end
% %     figure
% %     plot(y_conv)
%     y_bins_bands_uncom = fft(y_conv)/L;
%     y_bins_bands(runindex,:,1) = abs([y_bins_bands_uncom(1),y_bins_bands_uncom(2:L/2)+fliplr(y_bins_bands_uncom(L/2+2:L))]);
%     %[bestscore(runindex),bestbin(runindex)] = max(bins_bands(runindex,:,2));

end

%%
ave_num_pool = [1,2,4,8];
pd = zeros(1,length(ave_num_pool));
for ave_index = 1:length(ave_num_pool)
    clearvars bestscore bestbin
    ave_num = ave_num_pool(ave_index);
    if ave_num == 1
        for ii=1:runtimes
            [bestscore(ii),bestbin(ii)] = max(abs(bins_bands(ii,:,2)));
        end
    else
        for ii = 1:runtimes/ave_num
            [bestscore(ii),bestbin(ii)] = max(mean(bins_bands((ii-1)*ave_num+1:ii*ave_num,:,2)));
        end
    end

    pd(ave_index) = sum(bestbin ==3)/(runtimes/ave_num);
end
figure(2)
plot(ave_num_pool*Nn,pd,'x-');
hold on
grid on
xlabel('Number of Samples')
ylabel('Ident. Prob.')
ylim([0,1])

%% figure
figure(1)
for ii=1:2
    temp = squeeze(bins_bands(:,:,ii));
    ydata = mean(abs(temp),1);
    subplot(2,1,ii)
    plot(10*log10(ydata));
    grid on
end

% figure
figure(3)

temp = squeeze(y_bins_bands(:,:,1));
ydata = mean(abs(temp),1);
plot(10*log10(ydata));
grid on
