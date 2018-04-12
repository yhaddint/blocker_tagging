%  12/01/2015
%  test multi-resolution idea
%  Pseudo-spectrum is used. 

clear;clc;%clf;close all

blk_num = 31;
M = 1e2;
N = 32;

stat_num = 2e2;
P = M*N;
L = N*stat_num;

sig_pow_dB = [-3,-3,-30];
sig_pow = 10.^(sig_pow_dB/10);

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
%  other bands

Nn = 8*N;
L = N;
runtimes = 2;
bins_bands = zeros(runtimes,L,2);
for runindex = 1:runtimes
    phi = exp(1j*2*pi*rand(1,3));
    sig_sum = zeros(P,1);
    CAL = randi(2,P,3)*2-3;
    sig_cache = zeros(P,3);
    start_point = randi(P*M,1,3);
    for ii=1:3
        sig_cache(:,ii) = sig(start_point(ii):start_point(ii)+P-1,ii);
        if ii ~= 3
            sig_sum = sig_sum+2*real(sig_cache(:,ii).*phi(ii)).*CAL(:,1)*sqrt(sig_pow(ii));
        else
            sig_sum = sig_sum+2*real(sig_cache(:,ii).*phi(ii)).*CAL(:,2)*sqrt(sig_pow(ii));
        end
    end
    
    for ii=1:2
        sig_sum_corr = sig_sum.*CAL(:,ii);
        for ll=1:L
            sig_conv(ll) = sum(sig_sum_corr(1:Nn).*(sig_sum_corr((1:Nn)+ll)))/Nn;
        end
        bins_bands(runindex,:,ii) = abs(fft(sig_conv))/L;
    end
end
%% figure
for ii=1:2
    temp = squeeze(bins_bands(:,:,ii));
    ydata = mean(temp,1);
    subplot(2,1,ii)
    plot(10*log10([ydata(1),2*ydata(2:L/2)]));
    grid on
end